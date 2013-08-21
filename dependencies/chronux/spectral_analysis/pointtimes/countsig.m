function[H,P,M1,M2,N1,N2] = countsig(data1,data2,T1,T2,parametric,p,quiet)
% Give the program two spike data sets and one 
% or two time intervals and it will decide if  
% the counts are significantly different.      
% this is either with a non-parametric method  
% or with a sqrt transformation followed by a   
% t-test                                       
% Usage: [H,P,M1,M2,N1,N2] = countsig(data1,data2,T1,T2,parametric,p,quiet)
%                                              
% Input:                                       
% Note that all times have to be consistent. If data
% is in seconds, so must be sig and t. If data is in 
% samples, so must sig and t. The default is seconds.
%
% data1      - structure array of spike times (required)  
% data2      - structure array of spike times (required)  
% T1         - time interval (default all)     
% T2         - time interval (default T1)      
% parametric - 0 = non-parametric (Wilcoxon)   
%            - 1 = ttest on sqrt of counts     
%            - 2 = Poisson assumption          
%              (default = 0)                   
% p          - significance level (0.05)       
% quiet      - 1 = no display 0 = display      
%                                              
% Output:                                      
%                                              
% H          - 1 if different 0 if not         
% P          - prob of result if same          
% M1         - mean count for data1            
% M2         - mean count for data2            
% N1         - counts for data1                
% N2         - counts for data2                


if nargin < 2;error('I need 2 sets of spike data');end
data1=padNaN(data1); % create a zero padded data matrix from input structural array
data2=padNaN(data2); % create a zero padded data matrix from input structural array
data1=data1'; data2=data2'; % transpose data to get it into a form acceptable to Murray's routine
if nargin < 3, 
   T1 = [min(data1(:,1)) max(max(data1))]; 
end
if nargin < 4, 
   T2 = T1; 
end
if nargin < 5, 
   parametric = 0;
end  
if nargin < 6; p = 0.05;end
if nargin < 7; quiet = 0; end

if isempty(T1), 
   T1 = [min(data1(:,1)) max(max(data1))]; 
end
if isempty(T2) 
   T2 = T1; 
end
if isempty(parametric),
   parametric = 0;
end  
if isempty(p) 
   p = 0.05;
end
if isempty(quiet), 
    quiet = 0; 
end

NT1 = length(data1(:,1));
NT2 = length(data2(:,2));

if (NT1 < 4 || NT2 < 4) && parametric ~= 2,
  disp('Low number of trials : switch to Poisson test')
  parametric = 2;
end

if abs((T1(2)-T1(1)) - (T1(2)-T1(1))) > 10.^-6
  error('Time intervals for analysis are different')
end

% get counts...
N1=zeros(1,NT1);
for n=1:NT1
  N1(n) = length(find(data1(n,:) >=  T1(1) & ...
          data1(n,:) <=  T1(2) & ~isnan(data1(n,:))));
end
N2=zeros(1,NT2);
for n=1:NT2
  N2(n) = length(find(data2(n,:) >=  T2(1) & ...
          data2(n,:) <=  T2(2) & ~isnan(data1(n,:))));
end
M1 = mean(N1);
M2 = mean(N2);

% do non parametric test...

if parametric == 0
  [P H] = ranksum(N1,N2,p);
end

% parametric test (with stabilizing transform)...

%  use sqrt transformation from 
%  Cox and Lewis to make data more Gaussian
%  the statistical analysis of series of events pg 44

if parametric == 1
  X = sqrt(N1 +0.25);
  Y = sqrt(N2 +0.25);
  [H,P] = ttest2(X,Y,p,0);
end

%  Poisson test.  Use method from Zar
%  pg 580 Ed 3.  Z bahaves as a normal variate under
%  null of same process mean and Poisson processes

if parametric == 2
  X = sum(N1);
  Y = sum(N2);
  Z = abs(X-Y)./sqrt(X+Y);
  P = 2*(1-normcdf(Z));
  if P < p; H = 1;else H = 0;end
end  

if quiet == 0
  if H == 1
    disp('Counts are signifcantly different')
  else
    disp('Counts are not signifcantly different') 
  end
  disp(['Mean count for data1 = ' num2str(M1)])
  disp(['Mean count for data2 = ' num2str(M2)])
  disp(['P value = ' num2str(P)])
end











