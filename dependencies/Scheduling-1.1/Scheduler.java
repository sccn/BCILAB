import java.util.*;
import java.net.*;
import java.io.*;
import java.util.concurrent.*;
import java.util.concurrent.atomic.AtomicInteger;
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

public class Scheduler implements Runnable { 

	// create a new scheduler with a given set of workers (each specified as "host:port" string)
	public Scheduler(String[] workers) throws IOException { this(workers, "", 5, 1000, 10000); }

	// create a new scheduler with a given set of workers and a custom MATLAB scheduling policy
	//   all arguments must be passed; recommended defaults are 'par_reschedule_policy',3,3000,5000)
	public Scheduler(String[] workers, String custom_policy, int receiver_backlog, int receiver_timeout, int reschedule_interval) throws IOException {
		waiting = new TreeMap<Integer,String>();
		inflight = new ConcurrentHashMap<Integer,String>();
		finished = new ConcurrentHashMap<Integer,String>();
		events = new ConcurrentHashMap<Long,String>();

		exec = Executors.newCachedThreadPool();
		serv = new ServerSocket(0,receiver_backlog);
		serv.setSoTimeout(receiver_timeout);
		
		return_address = InetAddress.getLocalHost().getHostAddress() + ":" + Integer.toString(serv.getLocalPort());
		batchid = global_batchid.incrementAndGet();
		policy_interval = reschedule_interval;
		policy = custom_policy;
		shutdown = false;		

		// start threads (receiver, workers, policy)
		exec.execute(new ResultReceiver());
		for (String w : workers)
			exec.execute(new WorkerLink(w));
		exec.execute(this);
	}
	
	// submit the given tasks (string-formatted) for scheduling
	public synchronized void submit(String[] tasks) {
		for (String t : tasks)
			waiting.put(global_taskid.incrementAndGet(),t);
	}	

	// get results (string-formatted) that have been obtained since construction or last clear operation
	public synchronized String[] results() { return finished.values().toArray(new String[finished.values().size()]); }

	// test whether the scheduler is done with the tasks submitted since construction or last clear operation
	public synchronized boolean done() { return waiting.isEmpty() && inflight.isEmpty(); }

	// wait until the scheduler is done with the tasks submitted since the construction or last clear operation
	public synchronized void wait_until_done() {
		try {
			do {			
				wait();
			} while(!done());
		} catch (InterruptedException e) {}
	}

	// clear the scheduler's task & result sets 
	public synchronized void clear() { waiting.clear(); inflight.clear(); finished.clear(); events.clear(); batchid = global_batchid.incrementAndGet(); }

	// terminate background processing ASAP; renders the scheduler inoperational (only the results collected so far can be retrieved)
	public void terminate() { 
		shutdown = true; 
		try {
			serv.close();
		} catch (Exception e) {}
		exec.shutdownNow();
	}

	
	// --- internal data manipulation for waiting/inflight/result lists ---

	// atomically pick & assign a task to a worker (worker identified by a DataOutputStream)
	// returns id of the assigned task, or throws NoSuchElementException if the waiting queue is empty
	private Integer assign_task(Socket conn) throws NoSuchElementException, IOException {
		DataOutputStream out = new DataOutputStream(conn.getOutputStream());
		Integer id = null;
		String task = null;
		synchronized(this) {
			// get the next waiting entry (or throw if no such thing)
			id = waiting.firstKey();
			task = waiting.get(id);
			// and optimistically put it on the inflight queue
			inflight.put(id,task);
			waiting.remove(id);
		}
		try {
			// send it off
			out.writeInt(id);
			out.writeInt(return_address.length());
			out.writeBytes(return_address);    
			out.writeInt(task.length());
			out.writeBytes(task);
			out.flush();
			// check response
			DataInputStream in = new DataInputStream(conn.getInputStream());			
			int checksum = in.readInt();  
			if (checksum != id + return_address.length() + task.length())
				throw new IOException("Transmission Error.");
		} catch (IOException e) {
			// transmission error: put it back from inflight to waiting!
			synchronized(this) {
				waiting.put(id,task);
				inflight.remove(id);
			}
			throw e;
		}
		return id;
	}

	// atomically sign an incoming result off the inflight & waiting lists; 
	// notifies if both these lists become empty
	private synchronized void signoff_result(int id, String result) {
		// check that the id is actually valid (i.e. inflight or waiting (the latter if rescheduled in the meantime))
		boolean waits = waiting.containsKey(id);
		boolean flights = inflight.containsKey(id);
		if (flights || waits) {
			// valid: put into the results
			finished.put(id,result);
			// and remove from the lists
			if (waits)
				waiting.remove(id);
			if (flights)
				inflight.remove(id);
			if (waiting.isEmpty() && inflight.isEmpty())
				notifyAll();
		}
	}

	// atomically re-schedule a given in-flight task
	private synchronized void reschedule_task(int id) {
		if (inflight.containsKey(id)) {
			waiting.put(id,inflight.get(id));
			inflight.remove(id);
		}
	}

	
	// --- thread classes ---
	
	// thread that waits for incoming result-carrying connections, and spawns a ResultHandler thread for each of them
	private class ResultReceiver implements Runnable {  
		public void run() {    
			Socket conn = null;
			while (!shutdown) {
				try {
					// wait for a connection & hand off to result handler
					conn = serv.accept();
					exec.execute(new ResultHandler(conn));
					conn = null;
				} catch(InterruptedIOException e) {
					// accept timeout or tear-down requested: try to at least close the receiver connection
					try {
						if (conn != null)
							conn.close();
					} catch (Exception ex) { } 
				} catch(Exception e) {
					// register network or parsing fault
					if (!shutdown) {
						System.out.println("Exception while accepting results: " + e.getMessage());
						events.put(System.currentTimeMillis(),"error_accept");
					}
				}
			}
		}
	}

	// thread that receives a particular result
	private class ResultHandler implements Runnable {
		private Socket conn;
		public ResultHandler(Socket conn) { this.conn = conn; }
		// handle receiving and signing off of a result
		public void run() {
			int id = 0;
			try {
				// read the id & result
				DataInputStream in = new DataInputStream(conn.getInputStream());
				id = in.readInt();
				byte[] result = new byte[in.readInt()];
				in.readFully(result);	  
				// sign off the received data
				signoff_result(id,new String(result));	  
				// and register the event
				events.put(System.currentTimeMillis(),"received:" + Integer.toString(id) + ":" + conn.getInetAddress().getHostAddress());				
			} catch(IOException e) {
				// register network or parsing fault
				System.out.println("Exception while receiving results: " + e.getMessage());
				events.put(System.currentTimeMillis(),"error_receive:" + Integer.toString(id));
			} catch (Exception e) {
				System.out.println("Unexpected ResultHandler exception: " + e.getMessage());
			} finally {
				// clean up
				try {
					if (conn != null)
						conn.close();
				} catch (IOException e) {}
			}
		}
	}

	// thread that links to a given worker & assigns a task whenever the worker is ready to receive
	private class WorkerLink implements Runnable {
		private String address;  // address of our worker
		WorkerLink(String address) { this.address = address; }
		// handle the communication with the worker
		public void run() {
			try {
				String[] hostport = address.split(":");
				while (!shutdown) {
					Socket conn = null;
					try {
						if (!waiting.isEmpty()) {
							// ask for a slot in a worker 
							conn = new Socket(hostport[0],Integer.parseInt(hostport[1]));
							// acquired a slot: try to assign it a task 
							// (note: here we may find that all tasks have been scheduled to other workers in the meantime)
							try {
								Integer id = assign_task(conn);
								// register the successful assignment
								events.put(System.currentTimeMillis(),"assigned:" + Integer.toString(id) + ":" + conn.getInetAddress().getHostAddress());
							} catch (IOException e) {
								// checksum/transmission error...
								Thread.sleep(100);
							}
						} else {
							// queue empty: suspend for a while
							Thread.sleep(100);
						}
					} catch (NoSuchElementException e) {
						// queue empty; suspend for a while
						Thread.sleep(100);
					} catch (InterruptedException e) {
						// interrupted (tear-down requested)
						break;
					} catch (IOException e) {
						// random network problem: retry later, but do not spam worker with connect requests
						events.put(System.currentTimeMillis(),"error_assign:" + this.address);
						Thread.sleep(5000);
					} catch (Exception e) {
						// unexpected exception... possibly related to messed up data structures: wait
						System.out.println("Worker " + this.address + " encountered unexpected exception: "  + e.getMessage());
						Thread.sleep(5000);
					} finally {
						// clean up
						try {
							if (conn != null) 
								conn.close();
						} catch (IOException e) {}
					}
				}
			} catch (InterruptedException e) { }
		}
	}

	// scheduling policy thread
	public void run() {
		try {
			MatlabProxyFactory factory = new MatlabProxyFactory();
			MatlabProxy proxy = factory.getProxy();
			while (!shutdown) {
				try {
					// wait
					Thread.sleep(policy_interval);
					if (policy != null && !policy.isEmpty()) {
						// get a list of items to be re-scheduled (based on logged events & in-flight items)
						Object[] params = {(Integer)batchid,inflight.keySet().toArray(),waiting.keySet().toArray(),events.keySet().toArray(),events.values().toArray()};					
						Vector reschedule = (Vector)proxy.returningFeval(policy,1,params)[0];
						for (Object id : reschedule)
							reschedule_task((Integer)id);
					}
				} catch (InterruptedException e) {
					// interrupted (part of teardown procedure)
					break;
				} catch (Exception e) {
					System.out.println("Unexpected scheduling policy exception: " + e.getMessage());
				}
			}
			proxy.disconnect();
		} catch (MatlabConnectionException e) {
			System.out.println("Problem with the local MATLAB connection: " + e.getMessage());
		} catch (Exception e) {
			System.out.println("Unexpected scheduling policy thread exception: " + e.getMessage());
		}
	}


	private ExecutorService exec;  	// this one executes our threads
	private ServerSocket serv;  	// we wait for incoming results here
	private String return_address; 	// the target ("return") address of this scheduler, for the workers ("host:port")

	private SortedMap<Integer,String> waiting;  	// tasks (Id->StringRep) that are waiting to be scheduled; SYNCHRONIZED
	private ConcurrentMap<Integer,String> inflight; // tasks (Id->StringRep) that are (presumably) being worked on; SYNCHRONIZED
	private ConcurrentMap<Integer,String> finished; // the results (Id->StringRep) that have arrived so far
	private ConcurrentMap<Long,String> events; 	// a list of events (TimeMs->Event), for processing by the scheduling policy

	private boolean shutdown;	// shutdown notice for threads

	private String policy;	 	// re-scheduling policy function
	private int policy_interval;	// re-scheduling policy interval
	private int batchid;		// id of the current batch of operations
	private static AtomicInteger global_batchid = new AtomicInteger(0); // the global batch id (incremented)
	private static AtomicInteger global_taskid = new AtomicInteger(0); // the global task id (incremented)

	// --- test code ---

	public static void main(String[] args) throws IOException,InterruptedException {
		String[] tasks = {"sin(randn(10))","exp(randn(10))","sin(randn(10))","exp(randn(10))","sin(randn(10))","exp(randn(10))","sin(randn(10))","exp(randn(10))","sin(randn(10))","exp(randn(10))","sin(randn(10))","exp(randn(10))","sin(randn(10))","exp(randn(10))"};
		String[] workers = {"localhost:23547","localhost:23548"};
		Scheduler sched = new Scheduler(workers);
		sched.submit(tasks);
		sched.wait_until_done();
		String[] results = sched.results();
		for (String r : results)
			System.out.println(r);
		sched.terminate();
	}

}
