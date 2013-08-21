function varargout=pnet_remote(varargin)  
% PNET_REMOTE   - Evaluation of matlab expression in remote host PNET
%
%     Version: First includes in the tcp/udp/ip toolbox 2002-02-13
%              (C) 2002 Peter Rydesäter, GNU Public License 
%
%     This function uses PNET for nonblocking remote controll of other matlab 
%     session on this or other hosts. This function implements different client
%     and server calls in same function. By setting the server matlab session
%     in server mode with "PNET_REMOTE SERVER ..." its possible to connect to
%     it from other hosts and matlab sessions with "PNET_REMOTE con EVAL ..."
%     or exchange data with PUT or GET options. This remote controll package
%     uses its own (non standard) protocol over a TCP/IP connection.
%
%     Optional parameters are enclosed in [ ] in follwing syntax description.
% 
%  SYNTAX:
%  =======  
%
%  pnet_remote('server',[ port ])
%     
%     Starts to listen on specified port (is by default 5678) for a connection
%     and starts to serv for EVAL, PUT, GET... commands sent from remotehost.
%     It listens for a connection, servers it until it closes and then starts
%     to listens for a new connection to serv.
%   
%  pnet_remote(con,'serverat')
%  
%     Like the last but it starts to serv on an already estabished connection
%     until it is closed by the remote host. CON needs to be closed afterwards.
%
%  con=pnet_remote('connect',['host',[ port ] ])
%
%     Connects and returns a connection handler to a "PNET_REMOTE SERVER"
%     at specified HOST and PORT. Default value for host is 'localhost' and
%     for port 5678. A connection an also be established with
%     PNET('TCPCONNECT',...). 
%
%  pnet_remote(con,'close')
%
%     Sends a close command to remote host and closes the connection.
%     PNET(con,'CLOSE')  should also work even if it is not sending the close
%     command.
%  
%  pnet_remote(con,'close')
%
%     Exactly the same as PNET('CLOSEALL')
%  
%  pnet_remote(con,'eval','expression')
%
%     Non blocking evaluation of expression on remote host. The expression is
%     evaluated in callers namespace. The expression is a regular matlab
%     expression evaualted with EVALIN. This command is NONBLOCKING and
%     status of the evaluation can be detected with PNET_REMOTE(con,'STATUS')
%  
%  stat=pnet_remote(con,'status')
%
%     Returns evaluation status of a remote host/session. The returned value
%     is a string containing 'busy' 'error'  or 'ready'  depending on status
%     of evaluation. during evaluation 'busy' is returned, then 'error' or
%     'ready' is returened depending on if EVALIN succeded.
%
%  pnet_remote(con,'PUT','name1',expr1, 'name2',expr2.....)
%
%     Uploads variables to remote hosts evaulation workspace. Its specified
%     as a list of pairs of 'NAME' and EXPR where 'NAME' is the new name
%     of the uploaded variable in the remote workspace and EXPR is a local
%     variabel or expression.
%
%  [var1,var2,...]=pnet_remote(con,'GET','name1','name2'....)
%
%     Gets specified variables from remote workspace. Actually can 'name..'
%     be any expression that is remote evaluated _but_ this command blocks
%     until the result of the remotely evaluted expression/variable is  
%     recived.
%  
%  pnet_remote(con,'PUTSCRIPT','scriptname1','scriptname2',....)
%
%     Uploads local scripts to the remote host. The scripts are put on the 
%     remote hosts search path in a directory named PNET_PUTSCRIPTS.
%     'scriptname...' is the name used at calls. WHICH is used to detect
%     full path and extention of the script. If the remote hosts platform/OS is
%     different then the full extention on MEX files must be specified.
%
%  pnet_remote(con,'BREAK')
%  
%     Sends a break comand to the script running on the remote host. This break
%     comand can only be detected if the remote host repeatedly calls
%     PNET_REMOTE('GET_BREAK')  which will cause an error in the remote script.
%  
%  pnet_remote('GET_BREAK')  
%  
%     Used to be repeatedly called in scripts on the remote host. On a recived
%     BREAK or disconected connection an error will be generated that break
%     the running script.
%    
%  General alternative syntax(es):
%  
%       pnet_remote server port
%     
%     This is the same as the functional form PNET_REMOTE('SERVER',port) 
%     Numbers like connection handlers and port-numbers can be specified
%     as strings in most cases. This syntax should generaly work for all
%      variants of calls.
%
%       con_array=[con1 con2 con3 .....]
%       pnet_remote(con_array,.........)
%     
%     You can specify the conection handler as an array of connections and
%     same operation will be performed for all connections.
%     This syntax should generaly work for all variants of calls _IF_ they  
%     do not return anything.
%
 

%   This file(s) is part of the tcp_udp_ip toolbox (C) Peter Rydesäter et al.  
%   et al.  1998-2003 for running in MATLAB(R) as scripts and/or plug-ins.
%
%   This program is free software; you can redistribute it and/or modify
%   it under the terms of the GNU General Public License as published by
%   the Free Software Foundation; either version 2 of the License, or
%   (at your option) any later version.
%
%   This program is distributed in the hope that it will be useful,
%   but WITHOUT ANY WARRANTY; without even the implied warranty of
%   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%   GNU General Public License for more details.
%
%   You should have received a copy of the GNU General Public License
%   along with this program; if not, write to the Free Software
%   Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA
%
%   In addition, as a SPECIAL EXCEPTION, Peter Rydesäter, SWEDEN,
%   gives permission to link the code of this program with any library,
%   and distribute linked combinations of it. You must obey the GNU
%   General Public License in all respects for all of the code in the
%   tcp_udp_ip toolbox and any . If you modify any source file included,
%   you may extend this exception to your version of the file, but you are
%   not obligated to do so.  If you do not wish to do so, delete this exception
%   statement from your version. This exception makes it possible to use
%   pnet.c (.dll) as a plug-in as it is intended and let it be (dynamical)
%   linked to MATLAB(R) or any compiled stand alone application.
  
  persistent con;
  if isnumeric(varargin{1}),
    con=double(varargin{1});
    varargin=varargin(2:end);
    if length(con(:))>1,
      for con_n=con(:)',
	con_n
	pnet_remote(con_n,varargin{:});
      end
      return;
    end
  elseif ischar(varargin{1}),
    num=str2num(varargin{1});
    if length(num)==1,
      con=num(1);
      varargin=varargin(2:end);
    end
  end
  switch upper(varargin{1}),
   case 'CLOSEALL'
    pnet closeall;
    return;
   case 'CONNECT'
    rhost='127.0.0.1';
    rport=5678;
    con=-1;
    
    if length(varargin)>1, rhost=varargin{2}; end
    if length(varargin)>2, rport=varargin{3}; end
    while con==-1,
      con=pnet('tcpconnect',rhost,rport);
      if con==-1,
	disp(sprintf('CAN NOT CONNECT TO HOST: %s PORT: %d\nRETRY....',rhost,rport));
	pause(1);
      else
	disp(sprintf('CONNECTED TO HOST: %s PORT: %d !\n',rhost,rport));
      end
    end
    varargout{1}=con;
    return;
   case 'NEWSERVER'
    port=5678;
    if length(varargin)>1, port=varargin{2}; end
    if ~ischar(port), 
      port=sprintf('%d',port);
    end
    system(sprintf('echo "pnet_remote(''server'')" | matlab -nojvm -nosplash &',port));
    return;
   case 'GET_BREAK'
    global DEFAULT_CON__;
    if length(DEFAULT_CON__)==0,
      return;
    end
    if strcmp(local_status_str(DEFAULT_CON__),'break') | pnet(DEFAULT_CON__,'status')==0,
      pnet(DEFAULT_CON__,'printf','\n--error--\n');
      error 'Remote break';
    end
    return;
   case 'STATUS'
    varargout{1}=local_status_str(con);
    return;
   case 'WAITNOTBUSY'
    while strcmp(local_status_str(con),'busy'),
      pause(0.01);
    end
    return;
   case 'PUTSCRIPT'
    M=varargin(2:end);
    D={'PUTSCRIPT'};
    for n=1:length(M),
      try,
	f=[];file='';data='';file=which(M{n});f=fopen(file,'r');data=char(fread(f)),fclose(f);
      end
      if(length(file)>0 & length(data)>0),
	D{end+1}=file;D{end+1}=data;
      end
    end
    pnet(con,'printf','\n--remote--\n');
    local_status_str(con); % Flush status buffer. Keep last status in readbuffer
    pnet_putvar(con,D);
    return;
   case 'PUT'
    pnet(con,'printf','\n--remote--\n');
    local_status_str(con); % Flush status buffer. Keep last status in readbuffer
    pnet_putvar(con,varargin);
    return;
   case 'GET'
    pnet(con,'printf','\n--remote--\n');
    local_status_str(con); % Flush status buffer. Keep last status in readbuffer
    pnet_putvar(con,varargin);
    varargout=pnet_getvar(con);
    return;
   case 'EVAL'
    pnet_remote(con,'WAITNOTBUSY');
    pnet(con,'printf','\n--remote--\n');
    pnet_putvar(con,varargin);
    return;
   case 'CLOSE'
    pnet(con,'printf','\n--remote--\n');
    local_status_str(con); % Flush status buffer. Keep last status in readbuffer
    pnet_putvar(con,varargin);
    pnet(con,'close');
    return;
   case 'BREAK'
    pnet(con,'printf','\n--break--\n');
    return;
   case 'SERVER'
    port=5678;
    try, 
      port=varargin{2};
    end
    sock=pnet('tcpsocket',port);
    if 1,
      while 1,
	    disp(sprintf('WAIT FOR CONNECTION ON PORT: %d\n',port));
	    try,
	      con=[];
          con=pnet(sock,'tcplisten');
          con
	    catch,
	      disp 'Try:  "pnet closeall"  in all matlab sessions on this server.';
          disp ' ';
          error(lasterr);
	    end
	    try,
	      [rhost,rport]=pnet(con,'gethost');
	      disp(sprintf('START SERVING NEW CONNECTION FROM IP %d.%d.%d.%d port:%d',rhost,rport));
	      pnet_remote('SERVERNAMESPACE',con);
	    end
	    pnet(con,'close');
      end
    end
    pnet(sock,'close');
    disp('CLOSING SOCKET');
    return;
   case 'SERVERNAMESPACE'
    clear;
    pnet_remote(evalin('caller','con'),'serverat');
    return;    
   case 'SERVERAT'
    while pnet(con,'status'),
      okflag=1;
      if 1, %try,
	str='';
	drawnow;
	while strcmp(str,'--remote--')==0 & pnet(con,'status'),
	  str=pnet(con,'readline',1024);
	end
	if pnet(con,'status')==0, break; end
	C=pnet_getvar(con);
	pnet(con,'printf','\n--busy--\n');
	drawnow;
	switch upper(C{1}),
	 case 'EVAL'
	  global DEFAULT_CON__;
	  DEFAULT_CON__=con;
	  try
	    disp(['REMOTE EVAL>> ' C{2:min(2:end)}]);
	    evalin('caller',C{2:end},'okflag=0;');
	  catch
	    okflag=0;
	  end
	  DEFAULT_CON__=[];
	 case 'PUTSCRIPT'
	  C=C(2:end);
	  try, mkdir('PNET_PUTSCRIPTS'); end
	  try, addpath('PNET_PUTSCRIPTS'); end
	  for n=1:2:length(C),	  
	    disp(['REMOTE PUTSCRIPT>> ' C{n}]);
	    try,
	      [pa,na,ex]=fileparts(C{n});f=[];
	      f=fopen(['PNET_PUTSCRIPTS' filesep na ex],'w');fwrite(f,double(C{n+1}),'char');
	      fclose(f);
	    catch,
	      fclose(f);
	    end
	  end
	 case 'PUT'
	  C=C(2:end);
	  for n=1:2:length(C),	  
	    disp(['REMOTE PUT>> ' C{n}]);
	    try
	      assignin('caller',C{n:n+1});
	    catch
	      okflag=0;
	    end
	  end
	 case 'GET'
	  C=C(2:end);
	  R=cell(size(C));
	  for n=1:length(C),
	    disp(['REMOTE GET>> ' C{n}]);
	    try     
	      R{n}=evalin('caller',[C{n} ';']);
	    catch
	      okflag=0;
	    end
	  end
	  pnet_putvar(con,R);
	 case 'CLOSE'
	  return;
	end %END SWITCH
      end
      if okflag,
	pnet(con,'printf','\n--ready--\n');
      else
	pnet(con,'printf','\n--error--\n');
	disp(sprintf('\nERROR: %s\n',lasterr));
      end
    end
  end
  return;
 
% Function that returns and leavs last text line in buffer.
function stat=local_status_str(con)
  while 1, % Loop that finds, returns and leaves last text line in buffer.
    str=pnet(con,'read',    1024,'view','noblock');
    if length(find([str,' ']==char(10)))<=1,
      stat=pnet(con,'readline',1024,'view','noblock'); % The return 
      stat=stat(3:end-2);
      return;
    end
    dump=pnet(con,'readline',1024,'noblock'); % Then remove last line
  end
  return;
