function webserver_demo(port)
% WEBSERVER_DEMO - Demo WEBSERVER, returns a simple webpage to your browser
%
% Syntax:
%    webserver_demo
%  or
%    webserver_demo port
%
%  Version: 2002-02-01 for the tcpiptoolbox 2.x API
%
if(nargin==0),    port=8888; end
if(ischar(port)), port=str2num(port); end
sock=pnet('tcpsocket',port);
if(sock==-1), error('Specified TCP port is not possible to use now.'); end
pnet(sock,'setreadtimeout',1);
try,
  disp(sprintf(['Get a webpage with your browser at adress: http://localhost:%d\n' ...
	       'Or use proper hostname from an other computer.\n'],port));
  while 1,
    con=pnet(sock,'tcplisten');
    if( con~=-1 ),
      try,
	[ip,port]=pnet(con,'gethost');
	disp(sprintf('Connection from host:%d.%d.%d.%d port:%d\n',ip,port));
	pnet(con,'setreadtimeout',2);  % Avoid locking the server for slow/dead networks/browsers
	pnet(con,'setwritetimeout',1);  
	pnet(con,'readline');
	pnet(con,'printf','HTTP/1.1 200 OK\n');
	pnet(con,'printf','Content-Type: text/html\n\n');
	pnet(con,'printf','<html><HEAD><TITLE>WEBSERVER DEMO</TITLE></HEAD>\n<body><h1>WEBSERVER DEMO</h1>\n');
	str=sprintf('You are at host:%d.%d.%d.%d port:%d\n',ip,port);
	pnet(con,'printf','This webserver is a demo MATLAB(R) script using the tcpip-toolbox.<br><hr>%s<hr>\n',str);
	str=evalc('random_numbers=rand(5)');
	pnet(con,'printf','<h2>Some random numbers</h2>\n<pre>%s</pre>\n',str);
	pnet(con,'printf','<a href="./">Reload page with new numbers>></a>');
	pnet(con,'printf','<hr><I>(C) 2002 Peter Rydesäter, Mitthögskolan, Östersund, SWEDEN</I></body></html>\n');
      end
      pnet(con,'close');
      drawnow;
    end
  end
end
pnet(sock,'close');
