function gap = evalgap(fnc, fspec, dnorm, ww, uu, A, B, lambda)

[ff,gg]=evalloss(fnc,ww,uu,A,B);

fval = ff+lambda*sum(fspec(ww));
dval = evaldual(fnc,dnorm,-gg,A,B,lambda);

gap = (fval+dval)/fval;

function [fval,gg]=evalloss(fnc, ww, uu, A, B)

if ~isempty(uu)
  zz=A*ww+B*uu;
else
  zz=A*ww;
end

[fval, gg] =fnc.p(zz, fnc.args{:});


function dval = evaldual(fnc, dnorm, aa, A, B, lambda)

mm=length(aa);

if ~isempty(B)
  aa=aa-B*((B'*B)\(B'*aa));
end

[dnm,ishard] =dnorm(A'*aa);


if ishard && dnm>0
  aa  = min(1, lambda/dnm)*aa;
  dnm = 0; 
end

dval = fnc.d(aa, fnc.args{:})+dnm;


