function [x,N] = makesig(SigName,N)
% [x,N] = makesig(SigName,N) Creates artificial test signal identical to the
%     standard test signals proposed and used by D. Donoho and I. Johnstone
%     in WaveLab (- a matlab toolbox developed by Donoho et al. the statistics
%     department at Stanford University).
%
%    Input:  SigName - Name of the desired signal (Default 'all')
%                        'AllSig' (Returns a matrix with all the signals)
%                        'HeaviSine'
%                        'Bumps'
%                        'Blocks'
%                        'Doppler'
%                        'Ramp'
%                        'Cusp'
%                        'Sing'
%                        'HiSine'
%                        'LoSine'
%                        'LinChirp'
%                        'TwoChirp'
%                        'QuadChirp'
%                        'MishMash'
%                        'WernerSorrows' (Heisenberg)
%                        'Leopold' (Kronecker)
%            N       - Length in samples of the desired signal (Default 512)
%
%    Output: x   - vector/matrix of test signals
%            N   - length of signal returned
%
%    See also: 
%
%    References:
%            WaveLab can be accessed at
%            www_url: http://playfair.stanford.edu/~wavelab/
%            Also see various articles by D.L. Donoho et al. at
%            web_url: http://playfair.stanford.edu/
%
%Author: Jan Erik Odegard  <odegard@ece.rice.edu>
%This m-file is a copy of the  code provided with WaveLab
%customized to be consistent with RWT.

if(nargin < 1)
  SigName = 'AllSig';
  N = 512;
elseif(nargin == 1)
  N = 512;
end;
t = (1:N) ./N;
x = [];
y = [];
if(strcmp(SigName,'HeaviSine') | strcmp(SigName,'AllSig')),
  y = 4.*sin(4*pi.*t);
  y = y - sign(t - .3) - sign(.72 - t);
end;
x = [x;y];
y = [];
if(strcmp(SigName,'Bumps') | strcmp(SigName,'AllSig')),
  pos = [ .1 .13 .15 .23 .25 .40 .44 .65  .76 .78 .81];
  hgt = [ 4  5   3   4  5  4.2 2.1 4.3  3.1 5.1 4.2];
  wth = [.005 .005 .006 .01 .01 .03 .01 .01  .005 .008 .005];
  y = zeros(size(t));
  for j =1:length(pos)
    y = y + hgt(j)./( 1 + abs((t - pos(j))./wth(j))).^4;
  end 
end;
x = [x;y];
y = [];
if(strcmp(SigName,'Blocks') | strcmp(SigName,'AllSig')),
  pos = [ .1 .13 .15 .23 .25 .40 .44 .65  .76 .78 .81];
  hgt = [4 (-5) 3 (-4) 5 (-4.2) 2.1 4.3  (-3.1) 2.1 (-4.2)];
  y = zeros(size(t));
  for j=1:length(pos)
    y = y + (1 + sign(t-pos(j))).*(hgt(j)/2) ;
  end
end;
x = [x;y];
y = [];
if(strcmp(SigName,'Doppler') | strcmp(SigName,'AllSig')),
  y = sqrt(t.*(1-t)).*sin((2*pi*1.05) ./(t+.05));
end;
x = [x;y];
y = [];
if(strcmp(SigName,'Ramp') | strcmp(SigName,'AllSig')),
  y = t - (t >= .37);
end;
x = [x;y];
y = [];
if(strcmp(SigName,'Cusp') | strcmp(SigName,'AllSig')),
  y = sqrt(abs(t - .37));
end;
x = [x;y];
y = [];
if(strcmp(SigName,'Sing') | strcmp(SigName,'AllSig')),
  k = floor(N * .37);
  y = 1 ./abs(t - (k+.5)/N);
end;
x = [x;y];
y = [];
if(strcmp(SigName,'HiSine') | strcmp(SigName,'AllSig')),
  y = sin( pi * (N * .6902) .* t);
end;
x = [x;y];
y = [];
if(strcmp(SigName,'LoSine') | strcmp(SigName,'AllSig')),
  y = sin( pi * (N * .3333) .* t);
end;
x = [x;y];
y = [];
if(strcmp(SigName,'LinChirp') | strcmp(SigName,'AllSig')),
  y = sin(pi .* t .* ((N .* .125) .* t));
end;
x = [x;y];
y = [];
if(strcmp(SigName,'TwoChirp') | strcmp(SigName,'AllSig')),
  y = sin(pi .* t .* (N .* t)) + sin((pi/3) .* t .* (N .* t));
end;
x = [x;y];
y = [];
if(strcmp(SigName,'QuadChirp') | strcmp(SigName,'AllSig')),
  y = sin( (pi/3) .* t .* (N .* t.^2));
end;
x = [x;y];
y = [];
if(strcmp(SigName,'MishMash') | strcmp(SigName,'AllSig')),  
  % QuadChirp + LinChirp + HiSine
  y = sin( (pi/3) .* t .* (N .* t.^2)) ;
  y = y +  sin( pi * (N * .6902) .* t);
  y = y +  sin(pi .* t .* (N .* .125 .* t));
end;
x = [x;y];
y = [];
if(strcmp(SigName,'WernerSorrows') | strcmp(SigName,'AllSig')),
  y = sin( pi .* t .* (N/2 .* t.^2)) ;
  y = y +  sin( pi * (N * .6902) .* t);
  y = y +  sin(pi .* t .* (N .* t));
  pos = [ .1 .13 .15 .23 .25 .40 .44 .65  .76 .78 .81];
  hgt = [ 4  5   3   4  5  4.2 2.1 4.3  3.1 5.1 4.2];
  wth = [.005 .005 .006 .01 .01 .03 .01 .01  .005 .008 .005];
  for j =1:length(pos)
    y = y + hgt(j)./( 1 + abs((t - pos(j))./wth(j))).^4;
  end 
end;
x = [x;y];
y = [];
if(strcmp(SigName,'Leopold') | strcmp(SigName,'AllSig')),
  y = (t == floor(.37 * N)/N); 		% Kronecker
end;
x = [x;y];
y = [];

%  disp(sprintf('MakeSignal: I don*t recognize << %s>>',SigName))
%  disp('Allowable SigNames are:')
%  disp('AllSig'),
%  disp('HeaviSine'),
%  disp('Bumps'),
%  disp('Blocks'),
%  disp('Doppler'),
%  disp('Ramp'),
%  disp('Cusp'),
%  disp('Crease'),
%  disp('Sing'),
%  disp('HiSine'),
%  disp('LoSine'),
%  disp('LinChirp'),
%  disp('TwoChirp'),
%  disp('QuadChirp'),
%  disp('MishMash'),
%  disp('WernerSorrows'),
%  disp('Leopold'),
%end
