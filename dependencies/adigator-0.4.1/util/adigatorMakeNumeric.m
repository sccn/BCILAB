function y = adigatorMakeNumeric(x)
%function y = adigatorMakeNumeric(x)
% This routine turns a numeric object into an overloaded object.
% ----------------------- Input Information ----------------------------- %
% x:  numeric object
% --------------------- Output Information ------------------------------ %
% y:  overloaded object
%
% Copyright 2011-214 Matthew J. Weinstein and Anil V. Rao
% Distributed under the GNU General Public License version 3.0
global ADIGATOR
NUMvod = ADIGATOR.NVAROFDIFF;
VarID = ADIGATOR.VARINFO.COUNT;
[funcname,~] = cadafuncname();
func = struct('name',funcname,'size',size(x),'zerolocs',[],'value',x);
deriv = struct('name',cell(NUMvod,1),'nzlocs',cell(NUMvod,1));
y = cada(VarID,func,deriv);
ADIGATOR.VARINFO.LASTOCC(VarID,1) = VarID;
return
end