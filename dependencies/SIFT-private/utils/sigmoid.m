function v = sigmoid(params,range)

% Sigmoid creates a Sigmoid function using parameters in PARAMS and the 
% variable range.
% 
% V = SIGMOID(PARAMS,RANGE)
%
% PARAMS: a 3-vector, the entries of which are (in this order):
% amplitude, phase, slope, horizontal shift, vertical shift

amplitude = params(1);
Phase=params(2);
Slope=params(3);
hshift=params(4);
vshift=params(5);
%a=params(4);



v=vshift+amplitude*(1./(1+Phase*exp(-Slope*(range+hshift))));