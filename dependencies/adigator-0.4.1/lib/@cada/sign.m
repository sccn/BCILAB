function y = sign(x)
% CADA overloaded SIGN function: calls cadaunarymath
for Vcount = 1:length(x.deriv)
  if ~isempty(x.deriv(Vcount).nzlocs)
    warning('derivative of discontinuous sign function - making derivatives zero'); %#ok<WNTAG>
    x.deriv(Vcount).nzlocs = [];
    x.deriv(Vcount).name = [];
  end
end
y = cadaunarymath(x,1,'sign');