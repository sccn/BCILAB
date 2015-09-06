import java.util.*;
import java.net.*;
import java.io.*;
import java.util.concurrent.*;
import java.util.concurrent.atomic.AtomicInteger;
import java.text.DateFormat;
import java.text.SimpleDateFormat;
import matlabcontrol.*; // obtain the jar from http://code.google.com/p/matlabcontrol/ and add to the classpath


// Task scheduler for use with MATLAB (used by hlp_beginschedule.m / hlp_endschedule.m / hlp_schedule.m).
//
// * Constructed with a set of known workers and optional policy / control parameters.
// * the workers are assumed to be host:port addresses at which there may be processes listening
//   that receive (tagged) task messages, and forward (identically tagged) results to the collector
//   indicated in the task message.
// * the task message is formatted as (<>'s not part of the format):
//    <task_id><collectoraddress_length><collectoraddress><body_length><body>
//    where
//      <task_id> 					: tag/identifier of the task (int)
//      <collectoraddress_length>	: length, in bytes, of the data collector's address (int)
//      <collectoraddress>			: where to send the result for collection (string, formatted as host:port)
//      <body_length>				: length, in bytes, of the message body (int)
//      <body>						: task representation (string)
//  * the response message format emitted by the workers (and accepted by the Scheduler) is:
//    <task_id><body_length><body>
//
// The Scheduler ensures that the given tasks are assigned to available workers, and that results are collected;
// all processing happens in background threads, and is initiated via submit(), probed via done() / results(),
// and supports waiting for completion (of tasks submitted so far), via wait_until_done();
//
// Teardown of the Scheduler's background resources is initiated by executing terminate(), and it should be run
// once the Scheduler is no longer needed. If not invoked, the Scheduler's (spacious!) resources will remain in
// memory indefinitely, which leads to out-of-memory problems relatively quickly.
//
// Multiple schedulers can be accessing the same workers simultaneously, and the same scheduler can be fed with multiple
// submit()'s and reused for multiple independent scheduling rounds, by calling clear() to reset it into its initial state.
// Since a scheduler consumes (memory/network/cpu) resources, a teardown and later re-construction of a new scheduler should
// be considered between extended pauses; clear() does not reclaim these resources (and therefore allows more efficient restarts).
//
// A critical issue in every scheduler is the handling of tasks that were scheduled but lost in transmission (e.g. due to
// a worker crash). The chance of this happening rises with the number of workers involved. The re-scheduling of in-flight
// tasks is issued by the custom re-scheduling policy, which is a user-supplied MATLAB function. It receives an
// event log (tracking what was issued and/or received when, and what errors happened when) and a list of currently in-flight
// task ids, and emits a subset of these task ids for re-scheduling; it is invoked periodically by the scheduler.
//
// The scheduling policy needs a matlabcontrol jar file in the classpath; this functionality can be commented out without breaking
// the scheduler's main functionality, with the exception of re-scheduling tasks that were lost by workers or the network.
//
//                      Christian Kothe, Swartz Center for Computational Neuroscience, UCSD
//                      2010-08-28

class WorkerNotReadyException extends Exception {}

public class Scheduler implements Runnable {
    /** 
     * Create a new scheduler with a given set of workers and a custom MATLAB scheduling policy.
     * Note that all arguments must be passed.
     * @param workers: cell array of workers, in 'host:port' format
     * @param custom_policy: name of a MATLAB function that implements the rescheduling policy
     *                       (usually 'hlp_reschedule_policy')
     * @param result_acceptor: name of a MATLAB callback function that accepts results while the 
     * %                       scheduler is running (disabled if empty)
     * @param receiver_backlog: number of connections from workers returning results that can be 
     *                          queued up until they get processed (5 is a good choice)
     * @param receiver_timeout: When a worker that has previously announced that it is ready to 
     *                          transmit results does not start transmitting results within this 
     *                          time period (in ms) after pinged, the connection will be discarded 
     *                          and the worker has to re-connect (1000 is a good choice)
     * @param reschedule_interval: time between periodic checks if jobs have to be re-scheduled, in 
     *                             ms (5000 is a good choice)
     * @param verbosity_level: Verbosity level for the scheduling process (0-3)
     * @param size_hint: expected number of results that will be computed with this scheduler
     */
    public Scheduler(String[] workers, String custom_policy, String result_acceptor, int receiver_backlog, int receiver_timeout, int reschedule_interval, int verbosity_level, int size_hint) throws IOException {
        loglock = new Object();
        log("allocating data structures...");
        waiting = new TreeMap<Integer,byte[]>();
        inflight = new HashMap<Integer,byte[]>(workers.length*3);
        finished = new HashMap<Integer,byte[]>(size_hint);        
        events = new HashMap<Long,String>(size_hint*3);
        pass_out = new ConcurrentLinkedQueue<Integer>();
        rand = new Random();
        
        exec = Executors.newCachedThreadPool();
        serv = new ServerSocket(0,receiver_backlog);
        serv.setSoTimeout(receiver_timeout);
        
        return_address = InetAddress.getLocalHost().getHostAddress() + ":" + Integer.toString(serv.getLocalPort());
        batchid = global_batchid.incrementAndGet();
        policy_interval = reschedule_interval;
        policy = custom_policy;
        matlab_acceptor = result_acceptor;
        verbose = verbosity_level;
        shutdown = false;
        log("done.\n");
        
        try {
            log("initializing MATLAB proxy...\n");
            MatlabProxyFactory factory = new MatlabProxyFactory();
            proxy = factory.getProxy();
            log("MATLAB proxy is up.\n");
        } catch (MatlabConnectionException e) {
            errlog("MATLAB proxy could not be initialized: " + e.getMessage() + " (see matlabcontrol docs for troubleshooting).\n");
        }
        
        // start threads (receiver, workers, policy)
        log("launching threads...\n");
        exec.execute(new ResultReceiver(this));
        for (String w : workers)
            exec.execute(new WorkerLink(w,this));
        log("done launching threads.\n");
        exec.execute(this);        
    }

    // logging functions at various log levels
    public void errlog(String msg) { if (verbose>=0) writelog(msg); }
    public void log(String msg) { if (verbose>0) writelog(msg); }
    public void finelog(String msg) { if (verbose>1) writelog(msg); }
    public void veryfinelog(String msg) { if (verbose>2) writelog(msg); }
    public void writelog(String msg) { 
        DateFormat df = new SimpleDateFormat("MM/dd HH:mm:ss");
        Date now = Calendar.getInstance().getTime();
        synchronized(loglock) {
        System.out.print(df.format(now)+": " + msg);
        }
    }
    
    // submit the given tasks (string-formatted) for scheduling
    public synchronized void submit(Object[] tasks) {
        log("submit() -- adding " + Integer.toString(tasks.length) + " tasks to waiting list...\n");
        for (Object t : tasks) {
            Integer id = global_taskid.incrementAndGet();
            waiting.put(id,(byte[])t);
            log("  added id " + id.toString() + ".\n");
        }
        log("submit() completed.\n");
    }
    
    // get results that have been obtained since construction or last clear operation
    public synchronized byte[][] results() { return finished.values().toArray(new byte[finished.values().size()][]); }
    
    // test whether the scheduler is done with the tasks submitted since construction or last clear operation
    public synchronized boolean done() { return waiting.isEmpty() && inflight.isEmpty(); }
    
    // wait until the scheduler is done with the tasks submitted since the construction or last clear operation
    public synchronized void wait_until_done() {
        log("wait_until_done() -- now waiting until waiting and inflight are empty.\n");
        try {
            do {
                wait();
            } while(!done());
        } catch (InterruptedException e) {}
        log("wait_until_done() completed.\n");
    }
    
    // clear the scheduler's task & result sets
    public synchronized void clear() {
        log("clear() -- clearing data structures...");
        waiting.clear(); 
        inflight.clear(); 
        finished.clear(); 
        batchid = global_batchid.incrementAndGet();
        synchronized(events) { events.clear(); }
        log("clear() done.\n");
    }
    
    // terminate background processing ASAP; renders the scheduler inoperational (only the results collected so far can be retrieved)
    public void terminate() {
        log("terminate() -- initiating shutdown...\n");
        shutdown = true;
        clear();
        try {
            serv.close();
        } catch (Exception e) {}
        exec.shutdownNow();
        log("terminate() completed.\n");
        verbose = 0;
    }
    
    
    // --- internal data manipulation for waiting/inflight/result lists ---
    
    // atomically pick & assign a task to a worker (worker identified by a DataOutputStream)
    // returns id of the assigned task, or throws NoSuchElementException if the waiting queue is empty
    private Integer assign_task(Socket conn) throws NoSuchElementException, IOException, WorkerNotReadyException, Exception {
        log("    assign_tasks() -- attempting to assign a task to socket " + conn.getInetAddress().getHostAddress() + "...\n");
        log("      creating input/output streams for socket " + conn.getInetAddress().getHostAddress() + "...\n");
        DataInputStream in = new DataInputStream(conn.getInputStream());
        DataOutputStream out = new DataOutputStream(conn.getOutputStream());
        log("      done creating input/output streams for socket " + conn.getInetAddress().getHostAddress() + ".\n");
        Integer id = null;
        byte[] task = null;
        // check whether the worker is ready to accept a task
        try {
            log("      checking whether worker at socket " + conn.getInetAddress().getHostAddress() + " is ready...\n");
            int readycheck = in.readInt();
            if (readycheck != 12345)
                throw new WorkerNotReadyException();
            log("      worker at socket " + conn.getInetAddress().getHostAddress() + " confirmed as ready.\n");
        }
        catch (Exception e) {
            throw new WorkerNotReadyException();
        }
        // assign a waiting task id to the worker
        log("      attempting to obtain a task id for socket " + conn.getInetAddress().getHostAddress() + ".\n");
        synchronized(this) {
            // get the next waiting entry (or throw if no such thing)
            log("      looking for free task id for socket " + conn.getInetAddress().getHostAddress() + "...\n");
            id = waiting.firstKey();
            log("      found free task id (" + id.toString() + ") for socket " + conn.getInetAddress().getHostAddress() + "; moving from waiting to inflight list...\n");
            task = waiting.get(id);
            // and optimistically put it on the inflight queue
            inflight.put(id,task);
            waiting.remove(id);
            log("      id " + id.toString() + " moved to inflight list.\n");
        }
        try {
            log("      sending task id " + id.toString() + " over socket " + conn.getInetAddress().getHostAddress() + "...\n");
            // send it off
            out.writeInt(id);
            out.writeInt(return_address.length());
            out.writeBytes(return_address);
            out.writeInt(task.length);
            out.write(task);
            out.flush();
            log("      done sending task id " + id.toString() + " over socket " + conn.getInetAddress().getHostAddress() + ".\n");
            log("      confirming checksum for task id " + id.toString() + " over socket " + conn.getInetAddress().getHostAddress() + "...\n");
            // check response
            int checksum = in.readInt();
            if (checksum != id + return_address.length() + task.length)
                throw new IOException("Transmission Error.");
            log("      successfully confirmed checksum for task id " + id.toString() + " over socket " + conn.getInetAddress().getHostAddress() + ".\n");
        }
        catch (IOException e) {
            log("      got a transmission error trying to send task id " + id.toString() + " over socket " + conn.getInetAddress().getHostAddress() + ".\n");
            // transmission error: put it back from inflight to waiting!
            synchronized(this) {
                log("      putting task id " + id.toString() + " back to waiting list...\n");
                waiting.put(id,task);
                inflight.remove(id);
                log("      task id " + id.toString() + " has been put back on waiting list.\n");
            }
            throw e;
        }
        catch (Exception e) {
            // some other error
            log("      got an exception while trying to send task id " + id.toString() + " over socket " + conn.getInetAddress().getHostAddress() + ": " + e.getMessage() + ".\n");
        }
        return id;
    }
    
    // atomically sign an incoming result off the inflight & waiting lists;
    // notifies if both these lists become empty
    private void signoff_result(int id, byte[] result) {
        synchronized(this) {
            log("    now signing off result " + Integer.toString(id) + "...\n");
            // check that the id is actually valid (i.e. inflight or waiting (the latter if rescheduled in the meantime))
            boolean waits = waiting.containsKey(id);
            boolean flights = inflight.containsKey(id);
            if (flights || waits) {
                // id is valid...
                log("      result id " + Integer.toString(id) + " found in queue; adding to finished results...\n");
                // put into the results
                finished.put(id,result);
                log("      result with id " + Integer.toString(id) + " has been added to finished results.\n");
                if (!matlab_acceptor.isEmpty()) {
                    log("      result with id " + Integer.toString(id) + " being added to pass-out queue...\n");
                    pass_out.offer(id);
                    log("      result with id " + Integer.toString(id) + " successfully added to pass-out queue.\n");
                }     
                
                log("      result id " + Integer.toString(id) + " now being removed from waiting and inflight lists...\n");
                // and remove from the lists
                if (waits)
                    waiting.remove(id);
                if (flights)
                    inflight.remove(id);
                log("      result id "+ Integer.toString(id) + " has been removed from waiting and inflight lists...\n");
                if (waiting.isEmpty() && inflight.isEmpty()) {
                    log("waiting and inflight queues now empty.\n");
                    notifyAll();
                }
                log("    done signing off result " + Integer.toString(id) + ".\n");
            } else {
                log("      result id " + Integer.toString(id) + " is not currently waiting or inflight; ignoring.\n");
            }
        }
    }
    
    // atomically re-schedule a given in-flight task
    private synchronized void reschedule_task(int id) {
        log("    attempting to reschedule task " + Integer.toString(id) + ".\n");
        if (inflight.containsKey(id)) {
            log("      task " + Integer.toString(id) + " is currently in-flight; moving back to waiting...\n");
            waiting.put(id,inflight.get(id));
            inflight.remove(id);
            log("      task " + Integer.toString(id) + " has been rescheduled from in-flight to waiting.\n");
        } else {
            log("      task " + Integer.toString(id) + " could not be rescheduled from in-flight to waiting (reason: not currently in in-flight).\n");
        }
        log("    done rescheduling task " + Integer.toString(id) + ".\n");
    }
    
    
    // --- thread classes ---
    
    // thread that waits for incoming result-carrying connections, and spawns a ResultHandler thread for each of them
    private class ResultReceiver implements Runnable {
        private Runnable sched;
        public ResultReceiver(Runnable sched) { this.sched = sched; }
        public void run() {
            log("result receiver thread is initializing...\n");
            Socket conn = null;
            while (!shutdown) {
                try {
                    // wait for a connection & hand off to result handler
                    veryfinelog("result receiver is waiting for a connection...\n");
                    conn = serv.accept();
                    log("result receiver got a connection from " + conn.getInetAddress().getHostAddress() + "; issuing result handler...\n");
                    exec.execute(new ResultHandler(conn,this.sched));
                    conn = null;
                } catch(InterruptedIOException e) {
                    veryfinelog("result handler received interrupted IO exception: "+ e.getMessage() + "; waiting for retry...\n");
                    // accept timeout or tear-down requested: try to at least close the receiver connection
                    try {
                        if (conn != null) {
                            veryfinelog("result receiver closes connection...\n");
                            conn.close();
                            veryfinelog("result receiver has closed connection.\n");
                        }
                    } catch (Exception ex) { }
                } catch(Exception e) {
                    // register network or parsing fault
                    if (!shutdown) {
                        errlog("result receiver got an exception: " + e.getMessage() + "\n");
                        synchronized(events) { events.put(System.currentTimeMillis(),"error_accept"); }
                    }
                }
            }
            log("result receiver thread has shut down.\n");
        }
    }
    
    // thread that receives a particular result
    private class ResultHandler implements Runnable {
        private Socket conn;
        private String address;  // address of our worker
        private Runnable sched;
        public ResultHandler(Socket conn, Runnable sched) { this.conn = conn; this.sched = sched; this.address = conn.getInetAddress().getHostAddress(); }
        // handle receiving and signing off of a result
        public void run() {
            log("  result handler for " + this.address + " is initializing...\n");
            int id = 0;
            try {
                // read the id & result
                log("  result handler for " + this.address + " opens input stream...\n");
                DataInputStream in = new DataInputStream(conn.getInputStream());
                log("  result handler for " + this.address + " has opened input stream, now receiving data...\n");
                id = in.readInt();
                byte[] result = new byte[in.readInt()];
                in.readFully(result);
                log("  result for " + this.address + " (id=" + Integer.toString(id) + ") has been received.\n");
                log("  result handler for " + this.address + " now signing off result id " + Integer.toString(id) + " ...\n");
                // sign off the received data
                signoff_result(id,result);
                // and register the event
                synchronized(events) { events.put(System.currentTimeMillis(),"received:" + Integer.toString(id) + ":" + conn.getInetAddress().getHostAddress()); }
                log("  result id " + Integer.toString(id) + " is signed off.\n");
            } catch(IOException e) {
                // register network or parsing fault
                errlog("  IO error during result handling for " + this.address + " (" + e.getMessage() + ").\n");
                synchronized(events) { events.put(System.currentTimeMillis(),"error_receive:" + Integer.toString(id)); }
            } catch (Exception e) {
                log("  unexpected error during result handling for " + this.address + " (" + e.getMessage() + ").\n");
            } finally {
                // clean up
                log("  result handler is now shutting down...\n");
                try {
                    if (conn != null) {
                        veryfinelog("result handler closes connection...\n");
                        conn.close();
                        veryfinelog("result handler has closed connection.\n");
                    }
                } catch (IOException e) {}
            }
            log("  result handler has shut down...\n");
        }
    }
    
    // thread that links to a given worker & assigns a task whenever the worker is ready to receive
    private class WorkerLink implements Runnable {
        private String address;  // address of our worker
        private Runnable sched;
        WorkerLink(String address, Runnable sched) { this.address = address; this.sched = sched; }
        // handle the communication with the worker
        public void run() {
            try {
                log("  worker link to " + address + " initializing...\n");
                String[] hostport = address.split(":");
                while (!shutdown) {
                    Socket conn = null;
                    try {
                        boolean isempty;
                        synchronized(this.sched) { isempty = waiting.isEmpty(); }
                        if (!isempty) {
                            // ask for a slot in a worker
                            finelog("  worker link to " + address + " connecting...\n");
                            conn = new Socket(hostport[0],Integer.parseInt(hostport[1]));
                            log("  worker link to " + address + " (= " + conn.getInetAddress().getHostAddress() + ") established.\n");
                            // acquired a slot: try to assign it a task
                            // (note: here we may find that all tasks have been scheduled to other workers in the meantime)
                            try {
                                log("  worker link to " + address + " attempting to assign task...\n");
                                Integer id = assign_task(conn);
                                // register the successful assignment
                                synchronized(events) { events.put(System.currentTimeMillis(),"assigned:" + Integer.toString(id) + ":" + conn.getInetAddress().getHostAddress()); }
                                log("  worker link to " + address + " has successfully assigned task " + id.toString() + "\n");
                            } 
                            catch (WorkerNotReadyException e) {
                                finelog("  worker link to " + address + " found worker not ready; waiting for retry...\n");
                                // worker is not ready to accept...
                                Thread.sleep(500+rand.nextInt(1000));
                            }
                            catch (IOException e) {
                                log("  worker link to " + address + " encountered IO error (" + e.getMessage() + "); waiting for retry...\n");
                                // checksum/transmission error...
                                Thread.sleep(500+rand.nextInt(1000));
                            }
                        } else {
                            veryfinelog("  worker link to " + address + " found waiting-tasks queue empty; waiting to re-try...\n");
                            // queue empty: suspend for a while
                            Thread.sleep(500+rand.nextInt(1000));
                        }
                    } catch (NoSuchElementException e) {
                        // queue empty; suspend for a while
                        log("  worker link to " + address + " waiting for waiting-task queue to fill up again...\n");
                        Thread.sleep(500+rand.nextInt(1000));
                    } catch (InterruptedException e) {
                        // interrupted (tear-down requested)
                        break;
                    } catch (IOException e) {
                        // random network problem: retry later, but do not spam worker with connect requests
                        finelog("  worker link to " + address + " encountered network error (" + e.getMessage() + "); waiting for retry...\n");
                        synchronized(events) { events.put(System.currentTimeMillis(),"error_assign:" + this.address); }
                        Thread.sleep(2500+rand.nextInt(5000));
                    } catch (Exception e) {
                        // unexpected exception... possibly related to messed up data structures: wait
                        errlog("  worker link to " + address + " encountered unexpected error (" + e.getMessage() + "); waiting for retry...\n");
                        Thread.sleep(2500+rand.nextInt(5000));
                    } finally {
                        // clean up
                        try {
                            if (conn != null) {
                                veryfinelog("  worker link to " + address + " closing connection...\n");
                                conn.close();
                                veryfinelog("  worker link to " + address + " has closed connection.\n");
                            }
                        } catch (IOException e) {}
                    }
                }
            } catch (InterruptedException e) { }
            log("  worker link to " + address + " has shut down.\n");
        }
    }
    
    // housekeeping thread
    public void run() {
        try {
            log("entering housekeeping loop...\n");
            while (!shutdown) {
                try {
                    // wait
                    Thread.sleep(policy_interval);
                    if (policy != null && !policy.isEmpty()) {
                        // get a list of items to be re-scheduled (based on logged events & in-flight items)
                        Object[] params = new Object[5];
                        synchronized(this) {
                            params[0] = (Integer)batchid;
                            params[1] = inflight.keySet().toArray();
                            params[2] = waiting.keySet().toArray();
                        }
                        synchronized(events) {
                            params[3] = events.keySet().toArray();
                            params[4] = events.values().toArray();
                        }
                        Vector reschedule;
                        synchronized(proxy) { reschedule = (Vector)proxy.returningFeval(policy,1,params)[0]; }
                        for (Object id : reschedule)
                            reschedule_task((Integer)id);
                    }
                    // optionally process the pass-out queue by passing tasks out to the matlab_acceptor while storing only the tag in results
                    if (!matlab_acceptor.isEmpty()) {
                        while (true) {
                            Integer id = pass_out.poll();
                            if (id == null)
                                break;
                            // the tag acts as a reference to the actual result payload that we'll hand out to MATLAB
                            String tag = "tag__" + Integer.toString(rand.nextInt(2000000000));
                            // transfer tag+payload out to MATLAB
                            try {
                                log("      getting result content for task id " + Integer.toString(id) + ".\n");
                                byte[] result = null;
                                synchronized(this) { result = finished.get(id); }
                                log("      got result content for task id " + Integer.toString(id) + "; passing to MATLAB acceptor as tag " + tag + "...\n");
                                synchronized(proxy) { proxy.feval(matlab_acceptor,tag,result); }
                                log("      successfully passed result id " + Integer.toString(id) + " to MATLAB acceptor.\n");
                                // if (and only if) successful, we only store the tag in finished
                                log("      replacing result id " + Integer.toString(id) + " in finished by tag " + tag + "...\n");
                                synchronized(this) { finished.put(id, tag.getBytes()); }
                                log("      successfully replaced result id " + Integer.toString(id) + " in finished by tag " + tag + "...\n");
                            } catch (Exception e) {
                                errlog("failed to pass result id " + Integer.toString(id) + " out: " + e.getMessage() + "\n");
                            }
                        }
                    }
                } catch (InterruptedException e) {
                    // interrupted (part of teardown procedure)
                    break;
                } catch (Exception e) {
                    errlog("unexpected housekeeping thread exception: " + e.getMessage() + ".\n");
                }
            }
            log("exited housekeeping loop...\n");
            log("disconnecting from MATLAB proxy...\n");
            proxy.disconnect();
            log("MATLAB proxy is down.\n");
        } catch (Exception e) {
            errlog("Unexpected housekeeping thread exception: " + e.getMessage() + ".\n");
        }
    }
    

    private Random rand;                // random number generator for timings
    private MatlabProxy proxy;          // allows us to call MATLAB functions; do not call from within a block synchronized with Scheduler (can deadlock against MATLAB calling into e.g., done())
    private Object loglock;             // a lock to synchronize logging...

    private ExecutorService exec;  	// this one executes our threads
    private ServerSocket serv;  	// we wait for incoming results here
    private String return_address; 	// the target ("return") address of this scheduler, for the workers ("host:port")
    
    // these three lists are supposed to be only accessed from a synchronized(this) {} block or synchronized Scheduler function
    private SortedMap<Integer,byte[]> waiting;  // tasks (Id->TaskBytes) that are waiting to be scheduled
    private Map<Integer,byte[]> inflight; // tasks (Id->TaskBytes) that are (presumably) being worked on
    private Map<Integer,byte[]> finished; // the results (Id->TaskBytes) that have arrived so far

    private Map<Long,String> events;      // a list of events (TimeMs->Event), for processing by the scheduling policy
    private ConcurrentLinkedQueue<Integer> pass_out; // a list of task Id's in finished that shall be passed out to the matlab_acceptor at the next opportunity
    
    private boolean shutdown;	// shutdown notice for threads
    
    private String policy;          // re-scheduling policy function
    private String matlab_acceptor;	// result acceptor function
    private int policy_interval;	// re-scheduling policy interval
    private int verbose;		// verbosity level
    private int batchid;		// id of the current batch of operations
    private static AtomicInteger global_batchid = new AtomicInteger(0); // the global batch id (incremented)
    private static AtomicInteger global_taskid = new AtomicInteger(0); // the global task id (incremented)
    
    
    // --- test code ---
    
    public static void main(String[] args) throws IOException,InterruptedException {
        //System.out.print("main() -- running diagnostic tests...\n");
        //String[] tasks = {"sin(randn(10))","exp(randn(10))","sin(randn(10))","exp(randn(10))","sin(randn(10))","exp(randn(10))","sin(randn(10))","exp(randn(10))","sin(randn(10))","exp(randn(10))","sin(randn(10))","exp(randn(10))","sin(randn(10))","exp(randn(10))"};
        //String[] workers = {"localhost:23547","localhost:23548"};
        //Scheduler sched = new Scheduler(workers);
        //sched.submit(tasks);
        //sched.wait_until_done();
        //String[] results = sched.results();
        //for (String r : results)
        //	System.out.println(r);
        //sched.terminate();
        //System.out.print("main() -- done with diagnostic tests...\n");
    }
    
}
