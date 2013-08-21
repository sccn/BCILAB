% SPARSEBAYESDEMO  Simple demonstration of the SPARSEBAYES algorithm
%
%	SPARSEBAYESDEMO(LIKELIHOOD, DIMENSION, NOISETOSIGNAL)
%
% OUTPUT ARGUMENTS: None
% 
% INPUT ARGUMENTS:
% 
%	LIKELIHOOD		Text string, one of 'Gaussian' or 'Bernoulli'
%	DIMENSION		Integer, 1 or 2
%	NOISETOSIGNAL	An optional positive number to specify the
%					noise-to-signal (standard deviation) fraction.
%					(Optional: default value is 0.2).
% 
% EXAMPLES:
% 
%	SPARSEBAYESDEMO("Bernoulli",2)
%	SPARSEBAYESDEMO("Gaussian",1,0.5)
%
% NOTES: 
% 
% This program offers a simple demonstration of how to use the
% SPARSEBAYES (V2) Matlab software.
% 
% Synthetic data is generated from an underlying linear model based
% on a set of "Gaussian" basis functions, with the generator being
% "sparse" such that 10% of potential weights are non-zero. Data may be
% generated in an input space of one or two dimensions.
% 
% This generator is then used either as the basis for real-valued data with
% additive Gaussian noise (whose level may be varied), or for binary
% class-labelled data based on probabilities given by a sigmoid link
% function.
% 
% The SPARSEBAYES algorithm is then run on the data, and results and
% diagnostic information are graphed.
%

%
% Copyright 2009, Vector Anomaly Ltd
%
% This file is part of the SPARSEBAYES library for Matlab (V2.0).
%
% SPARSEBAYES is free software; you can redistribute it and/or modify it
% under the terms of the GNU General Public License as published by the Free
% Software Foundation; either version 2 of the License, or (at your option)
% any later version.
%
% SPARSEBAYES is distributed in the hope that it will be useful, but WITHOUT
% ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
% FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for
% more details.
%
% You should have received a copy of the GNU General Public License along
% with SPARSEBAYES in the accompanying file "licence.txt"; if not, write to
% the Free Software Foundation, Inc., 51 Franklin St, Fifth Floor, Boston,
% MA 02110-1301 USA
%
% Contact the author: m a i l [at] m i k e t i p p i n g . c o m
%
function SparseBayesDemo(likelihood_, dimension, noiseToSignal)

% Fix the random seed for reproducibility of results
% 
rseed	= 1;
rand('state',rseed);
randn('state',rseed);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% --- VALIDATE INPUT ARGUMENTS ---
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% (1) likelihood
% 
LIKELIHOOD	= SB2_Likelihoods(likelihood_);
%
% (2) dimension
%
if dimension~=1 && dimension~=2
  error('Specified dimension should be 1 or 2')
end
%
% Set up default for "noiseToSignal" variable.
% For ease of use, we'll just ignore it in the case of a non-Gaussian
% likelihood model.
%
if ~exist('noiseToSignal','var')
  noiseToSignal	= 0.2;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% --- SET UP DEMO PARAMETERS ---
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% Experiment with these values to vary the demo if you wish
%
if dimension==1
  N	= 100;	% Number of points
else
  N	= 900;	% Gives a nice square grid of decent size
end
%
basisWidth	= 0.05;		% NB: data is in [0,1]
%
% Define probability of a basis function NOT being used by the generative
% model. i.e. if pSparse=0.90, only 10% of basis functions (on average) will
% be used to synthesise the data.
% 
pSparse		= 0.90;
iterations	= 500;
%
% Heuristically adjust basis width to account for 
% distance scaling with dimension.
% 
basisWidth	= basisWidth^(1/dimension);
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% --- SYNTHETIC DATA GENERATION ---
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% First define the input data over a regular grid
% 
if dimension==1
  X	= [0:N-1]'/N;
else
  % dimension is 2
  sqrtN		= floor(sqrt(N));
  N			= sqrtN*sqrtN;
  x			= [0:sqrtN-1]'/sqrtN;
  [gx gy]	= meshgrid(x);
  X			= [gx(:) gy(:)];
end
%
% Now define the basis 
% 
% Locate basis functions at data points
% 
C	= X;
%
% Compute ("Gaussian") basis (design) matrix
% 
BASIS	= exp(-distSquared(X,C)/(basisWidth^2));
%
%
% Randomise some weights, then make each weight sparse with probability
% pSparse 
% 
M			= size(BASIS,2);
w			= randn(M,1)*100 / (M*(1-pSparse));
sparse		= rand(M,1)<pSparse;
w(sparse)	= 0;
%
% Now we have the basis and weights, compute linear model
% 
z			= BASIS*w;
%
% Finally generate the data according to the likelihood model
% 
switch (LIKELIHOOD.InUse)
 case LIKELIHOOD.Gaussian,
  % Generate our data by adding some noise on to the generative function
  noise		= std(z) * noiseToSignal;
  Outputs	= z + noise*randn(N,1);
  %
 case LIKELIHOOD.Bernoulli,
  % Generate random [0,1] labels given by the log-odds 'z'
  Outputs	= double(rand(N,1)<SB2_Sigmoid(z));
end
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% --- SET UP GRAPHING PARAMETERS ---
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
fRows		= 2;
fCols		= 3;
%
SP_DATA		= 1;
SP_LIKELY	= 2;
SP_LINEAR	= 3;
SP_COMPARE	= 4;
SP_WEIGHTS	= 5;
SP_GAMMA	= 6;
%
figure(1)
clf
whitebg(1,'k')
TITLE_SIZE	= 12;
%
% Plot the data 
% 
subplot(fRows,fCols,SP_DATA)
if dimension==1
  plot(X,Outputs,'w.')
else
  plot3(X(:,1),X(:,2),Outputs,'w.')
end
t_	= sprintf('Generated data (%d points)', N);
title(t_,'FontSize',TITLE_SIZE)
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% --- SPARSE BAYES INFERENCE SECTION ---
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% The section of code below is the main section required to run the
% SPARSEBAYES algorithm.
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Set up the options:
% 
% - we set the diagnostics level to 2 (reasonable)
% - we will monitor the progress every 10 iterations
% 
OPTIONS		= SB2_UserOptions('iterations',iterations,...
							  'diagnosticLevel', 2,...
							  'monitor', 10);
%
% Set initial parameter values:
% 
% - this specification of the initial noise standard deviation is not
% necessary, but included here for illustration. If omitted, SPARSEBAYES
% will call SB2_PARAMETERSETTINGS itself to obtain an appropriate default
% for the noise (and other SETTINGS fields).
% 
SETTINGS	= SB2_ParameterSettings('NoiseStd',0.1);
%
% Now run the main SPARSEBAYES function
%
[PARAMETER, HYPERPARAMETER, DIAGNOSTIC] = ...
    SparseBayes(likelihood_, BASIS, Outputs, OPTIONS, SETTINGS)
%
% Manipulate the returned weights for convenience later
%
w_infer						= zeros(M,1);
w_infer(PARAMETER.Relevant)	= PARAMETER.Value;
%
% Compute the inferred prediction function
% 
y				= BASIS*w_infer;
%
% Convert the output according to the likelihood (i.e. apply link function)
% 
switch LIKELIHOOD.InUse
 case LIKELIHOOD.Gaussian
  y_l	= y;
 case LIKELIHOOD.Bernoulli
  y_l	= double(SB2_Sigmoid(y)>0.5);
end
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% --- PLOT THE RESULTS ---
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Likelihood trace (and Gaussian noise info)
% 
subplot(fRows,fCols,SP_LIKELY)
lsteps	= length(DIAGNOSTIC.Likelihood);
plot(1:lsteps, DIAGNOSTIC.Likelihood, 'g-')
set(gca,'Xlim',[0 lsteps+1])
title('Log marginal likelihood trace','FontSize',TITLE_SIZE)
if LIKELIHOOD.InUse==LIKELIHOOD.Gaussian
  ax	= axis;
  dx	= ax(2)-ax(1);
  dy	= ax(4)-ax(3);
  t_	= sprintf('Actual noise:   %.5f', noise);
  text(ax(1)+0.1*dx,ax(3)+0.6*dy,t_,'FontName','Courier')
  t_	= sprintf('Inferred noise: %.5f', 1/sqrt(HYPERPARAMETER.beta));
  text(ax(1)+0.1*dx,ax(3)+0.5*dy,t_,'FontName','Courier')
end
%
% Compare the generative and predictive linear models
% 
subplot(fRows,fCols,SP_LINEAR)
if dimension==1
  plot(X,z,'w-','linewidth',4);
  hold on
  plot(X,y,'r-','linewidth',3);
  hold off
else
  mesh(gx,gy,reshape(z,size(gx)),'edgecolor','w','facecolor','w')
  hold on
  mesh(gx,gy,reshape(y,size(gx)),'edgecolor','r','facecolor','r')
  hold off
  light
end
title('Generative function and linear model','FontSize',TITLE_SIZE)
legend('Actual','Model','Location','NorthWest')
legend('boxoff')
%
% Compare the data and the predictive model (post link-function)
% 
subplot(fRows,fCols,SP_COMPARE)
if dimension==1
  plot(X,Outputs,'w.','linewidth',4);
  hold on
  plot(X,y_l,'r-','linewidth',3);
  plot(X(PARAMETER.Relevant),Outputs(PARAMETER.Relevant),'yo')
  hold off
else
  plot3(X(:,1),X(:,2),Outputs,'w.')
  hold on
  mesh(gx,gy,reshape(y_l,size(gx)),'edgecolor','r','facecolor','r')
  hold off
  light
end
title('Data and predictor','FontSize',TITLE_SIZE)
%
% Show the inferred weights
% 
subplot(fRows,fCols,SP_WEIGHTS)
h	= stem(w_infer,'filled');
set(h,'Markersize',3,'Color','r')
set(gca,'Xlim',[0 N+1])
t_	= sprintf('Inferred weights (%d)', length(PARAMETER.Relevant));
title(t_,'FontSize',TITLE_SIZE)
%
% Show the "well-determinedness" factors
% 
subplot(fRows,fCols,SP_GAMMA)
bar(DIAGNOSTIC.Gamma,'g')
axis([0 length(PARAMETER.Relevant)+1 0 1.1])
title('Well-determinedness (gamma)','FontSize',TITLE_SIZE)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% Support function to compute basis
%
function D2 = distSquared(X,Y)
%
nx	= size(X,1);
ny	= size(Y,1);
%
D2 = (sum((X.^2), 2) * ones(1,ny)) + (ones(nx, 1) * sum((Y.^2),2)') - ...
     2*X*Y';

