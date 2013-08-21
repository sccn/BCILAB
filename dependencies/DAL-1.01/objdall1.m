% objdall1 - objective function of DAL with L1 regularization
%
% Copyright(c) 2009 Ryota Tomioka
% This software is distributed under the MIT license. See license.txt

function varargout=objdall1(aa, info, prob, ww, uu, A, AT, B, lambda, eta)

m = length(aa);
n = length(ww);
vv = AT(aa)+ww/eta;

Ip = find(vv>lambda);
In = find(vv<-lambda);

if nargout<=3
  [floss, gloss, hmin]=prob.floss.d(aa, prob.floss.args{:});
else
  [floss, gloss, hloss, hmin]=prob.floss.d(aa, prob.floss.args{:});
end


vsth = l1_softth(vv,lambda);


fval = floss+0.5*eta*sum(vsth.^2);
if ~isempty(uu)
  u1   = uu/eta+B'*aa;
  fval = fval + 0.5*eta*sum(u1.^2);
end

varargout{1}=fval;


if nargout<=2
  varargout{2} = info;
else
  gg  = gloss+eta*(A(vsth));
  soc = sum((vsth-ww/eta).^2);
  if ~isempty(uu)
    gg  = gg+eta*(B*u1);
    soc = soc+sum((B'*aa).^2);
  end

  if soc>0
    info.ginfo = norm(gg)/(sqrt(eta*hmin*soc));
  else
    info.ginfo = inf;
  end
  varargout{2} = gg;

  if nargout==3
    varargout{3} = info;
  else
    I = sort([Ip; In]);
    len = length(I);
    F = sparse(I,1:len, ones(1,len),n,len);
    AF = A(F);

    switch(info.solver)
     case 'cg'
      prec=hloss+spdiag(eta*sum(AF.^2,2));
      if ~isempty(uu)
        prec =prec+spdiag(eta*sum(B.^2,2));
      end
      varargout{3} = struct('hloss',hloss,'AF',AF,'I',I,'n',n,'prec',prec,'B',B);
     otherwise
      if length(I)>0
        varargout{3} = hloss+eta*AF*AF';
% $$$         sp=svd(varargout{3});
% $$$         cond=max(sp)/min(sp);
% $$$         fprintf('cond=%g\n',cond);
      else
        varargout{3} = hloss;
      end
      if ~isempty(uu)
        varargout{3} = varargout{3}+eta*B*B';
      end
    end
    varargout{4} = info;
  end
end

