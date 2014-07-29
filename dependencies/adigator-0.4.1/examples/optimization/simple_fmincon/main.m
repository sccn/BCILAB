% In this example we use MATLAB's fmincon function to solve the objective
% function in the file objfun.m with the inequality constraints in function
% confun.m. If the user does not have MATLAB's optimization package
% installed, then fmincon will not be called, but the derivative file
% creation will still be performed.
% Copyright 2011-2014 Matthew J. Weinstein and Anil V. Rao
% Distributed under the GNU General Public License version 3.0
clear all
xi = [-1,1];
solveflag = exist('fmincon','file');

% ------------------- Solve Without Derivatives ------------------------- %
if solveflag
  tic;
  options = optimset('Algorithm','interior-point','Display','iter');
  [x0,fval0] = ...
    fmincon(@objfun,xi,[],[],[],[],[],[],@confun,options);
  time0 = toc;
end

% ------------------- Solve With 1st Derivatives ------------------------ %
% Need to take derivatives of objfun and confun
gx = adigatorCreateDerivInput([1 2],'x'); % Deriv Input for x
adigator('objfun',{gx},'objfun_x',adigatorOptions('overwrite',1)); % Create Obj Deriv File
adigator('confun',{gx},'confun_x',adigatorOptions('overwrite',1)); % Create Con Deriv File
% We also need to make wrapper functions for each of these - these are made
% in the files objfungrad and confunjac.
if solveflag
  tic;
  options = optimset('Algorithm','interior-point');
  options = optimset(options,'GradObj','on','GradConstr','on','Display','iter');
  [x1,fval1] = fmincon(@objfungrad,xi,[],[],[],[],[],[],...
    @confungrad,options);
  time1 = toc;
end

% -------------- Solve With 2nd Derivatives - Option 1 ------------------ %
% To solve with 2nd derivatives we need to create a file which computes the
% Langrangian Hessian - The first way which we can do this is to simply
% build the Langrangian and then take the derivative of it twice - we build
% this function in the file lagrangian1.
glambda = adigatorCreateAuxInput([2 1]);
adigator('lagrangian1',{gx,glambda},'lagrangian1_x',adigatorOptions('overwrite',1)); % First deriv of lagrangian1
% For second deriv we need to make the input a structure since thats the
% input for our first deriv file
out = adigator('lagrangian1_x',{struct('f',gx,'dx',[1;1]),glambda},'lagrangian1_xx',...
  adigatorOptions('overwrite',1,'comments',0)); 
% We now need to make a wrapper for this, we called it laghess1.m
if solveflag
  tic;
  options = optimset('Algorithm','interior-point',...
    'Display','iter','GradObj','on','GradConstr','on',...
    'Hessian','user-supplied','HessFcn',@laghess1);
  [x2,fval2] = fmincon(@objfungrad,xi,[],[],[],[],[],[],...
    @confungrad,options);
  time2 = toc;
end

% -------------------- Solve With 2nd Derivatives - Option 2 ------------ %
% Rather than doing what we did in the last section, we can note that we've
% already created files for the objective gradient and constraint jacobian,
% thus we could create a file which builds the Lagrangian gradient, where
% GL = GO + lambda.'*JC, where 
% GL = lagrangian gradient
% GO = objective gradient  - computed in file objfun_x
% JC = constraint jacobian - computed in file confun_x
% So we make a file called laggrad and then differentiate this file
adigator('laggrad',{gx,glambda},'laggrad_x',adigatorOptions('overwrite',1));
% We now need to make a wrapper for this so that we can build the Hessian -
% we call this file laghess2
if solveflag
  tic;
  options = optimset('Algorithm','interior-point',...
    'Display','iter','GradObj','on','GradConstr','on',...
    'Hessian','user-supplied','HessFcn',@laghess2);
  [x3,fval3] = fmincon(@objfungrad,xi,[],[],[],[],[],[],...
    @confungrad,options);
  time3 = toc;
end

% Display Solve Times
if solveflag
  fprintf(['Solve time using no derivatives',...
    ': ',num2str(time0),'\n'])
  fprintf(['Solve time using 1st derivatives',...
    ': ',num2str(time1),'\n'])
  fprintf(['Solve time using 2nd derivatives',...
    ' (option1): ',num2str(time2),'\n'])
  fprintf(['Solve time using 2nd derivatives',...
    ' (option2): ',num2str(time3),'\n'])
end