function [fval_true,TVval,FIDval] = funcval(U,aTV,aL1,B,picks)

global Ux Uy PsiTU

[m,n] = size(U);

fval_true = 0;

if aL1 > 0
    fval_true = fval_true + aL1* sum(abs(PsiTU(:)));
end

if isempty(Ux) || isempty(Uy)
   [m,n] = size(U);
   Ux = [diff(U,1,1);U(1,:)-U(end,:)];
   Uy = [diff(U,1,2),U(:,1)-U(:,end)];
end


if aTV > 0
    temp = sum(sum(sqrt(Ux.^2 + Uy.^2)));
    fval_true = fval_true + aTV*temp;
    if nargout > 1
        TVval = temp;
    end
end

fu = fft2(U);
temp = fu(picks)/sqrt(m*n) - B;
temp = norm([real(temp);imag(temp)]);
temp = .5*temp^2;
fval_true = fval_true + temp;%/(m*n);
if nargout > 2
    FIDval = temp;
end


