function[N,B,E] = isi(data,T,err,Nbins,plt)
% Calculate the inter-spike-interval histogram                 
%    Usage: [N,B,E] = isi(data,T,err,Nbins,plt)
%    
% Input:                                                       
% Note that all times have to be consistent. 
%
% data   - structure array of spike times  (required)             
% T      - time interval of interest (default all)             
% err    - 0 for no error bars, 1 for jackknife errors
% 
% Nbins  - number of bins in the isi                           
%                                                              
% Output:                                                      
%                                                              
% N      - count in bins                                       
% B      - bin centres                                            
% E      - errorbar (this is 2 sig deviation                   
%          calculated using a jackknife over trials)           


if nargin < 1; error('I need data!'); end
data=padNaN(data); % create a zero padded data matrix from input structural array
data=data'; % transposes data to get it in a form compatible with Murray's routine
if nargin < 2; T = [min(data(:,1)) max(max(data))]; end
if nargin < 3; err = 0;end
if nargin < 4; Nbins = -1; end
if nargin < 5; plt = 'r'; end

if isempty(T); T = [min(min(data)) max(max(data))]; end
if isempty(err); err = 0;end
if isempty(Nbins); Nbins = -1; end
if isempty(plt); plt = 'r'; end

%  get the number of intervals in each trial and the indices of spike times
%  that are kept

NT = length(data(1,:)); % number of trials
NI=zeros(1,NT);
index(1:NT)=struct('keep',[]);
for n=1:NT
  indx = find(data(:,n) >=  T(1) & data(:,n) <=  T(2) ... 
                                 & ~isnan(data(:,n)));
  if isempty(indx)
    NI(n) = 0;
  else
    NI(n) = length(indx)-1;
    index(n).keep=indx;
  end 
end


% calculate intervals...

I = zeros(NT,max(NI));
IT = [];
for n=1:NT
  I(n,1:NI(n)) = diff(data(index(n).keep,n));
  IT = [IT I(n,1:NI(n))];
end

Mx = max(IT);
if Nbins == -1
  Nbins = floor(sum(NI)/30);
  Med = median(IT);
  Nbins = max(floor(Nbins*Mx/Med),10);
end

B = linspace(0,Mx,Nbins);

N = zeros(NT,Nbins);
for n=1:NT
  N(n,:) = hist(I(n,1:NI(n)),B);
end

% answer...

if NT > 1;Ns = sum(N)/NT;else Ns = N;end
if ~strcmp(plt,'n')
  bar(B,NT*Ns);
end

% Jackknife iver trials to estimate std...

if NT > 4 && err == 1
  MN = 0;
  SN = 0;
  for n=1:NT
    JK = (NT*Ns - N(n,:))/(NT-1);
    MN = MN + JK;
    SN = SN + JK.^2;   
  end  
  MN = MN/NT;
  SN = SN/NT;
  E = sqrt((NT-1)*(SN - MN.^2));
  if ~strcmp(plt,'n')
    hold on
    errorbar(B,NT*Ns,NT*2*E,'r-')
    hold off
  end
end
N = NT*Ns;








