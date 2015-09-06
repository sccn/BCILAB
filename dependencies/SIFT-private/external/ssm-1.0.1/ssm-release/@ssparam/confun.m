function fun = confun(param)

%@SSPARAM/CONFUN 

% (c) 2006-2007 Jyh-Ying Peng ´^´¼·­
% $Revision 1.0.0 $ $Generated: 2007/09/04 $

%%%%%%%% TODO: Use constraint function?
N   = size(param.value, 1);

if N > 1, fun = @bounceN; else fun = @bounce; end

n = length(param.group)-1;
t = zeros(1, n);
for i = 1 : n
    switch param.transform{i}
        case 'covariance', t(i) = 1;
        case 'ar>=3', t(i) = 2;
        case 'ma>=3', t(i) = 3;
    end
end
g = param.group;

%     function [C Ceq] = con(XX)
%         C   = []; % C <= 0
%         Ceq = 0; % Ceq == 0
%         for i = 1 : n
%             if t(i) == 1 % covariance
%                 X = XX(g(i)+1:g(i+1));
%                 p = floor(realsqrt(2*length(X)));
%                 Q = diag(exp(2*X(1:p)));
%                 Q(logical(tril(ones(p))-eye(p))) = X(p+1:end)./sqrt(1+X(p+1:end).^2);
%                 for j = 2:p, for k = 1:j-1, Q(j, k) = sqrt(Q(j, j)*Q(k, k))*Q(j, k); end, end
%                 Q = Q + Q' - diag(diag(Q));
%                 C = [C; -det(Q)];
%             elseif t(i) == 2 % ar>=3
%                 C = [C; abs(roots([1 -XX(g(i)+1:g(i+1))]))-1];
%             elseif t(i) == 3 % ma>=3
%                 C = [C; abs(roots([1 XX(g(i)+1:g(i+1))]))-1];
%             end
%         end
%     end

    function XX = bounce(XX)
        for i = 1 : n
            switch t(i)
                case 1 % covariance
                    %%%% TODO: bounce negative covariance?
                    % warning('ssm:ssparam:NegCov', 'Negative covariance.\n');
                case 2 % ar>=3
                    r = roots([1 -XX(g(i)+1:g(i+1))]);
                    v = abs(r) > 1;
                    if any(v)
                        r(v) = 1./r(v);
                        X = poly(r);
                        XX(g(i)+1:g(i+1)) = -X(2:end);
                    end
                case 3 % ma>=3
                    r = roots([1 XX(g(i)+1:g(i+1))]);
                    v = abs(r) > 1;
                    if any(v)
                        r(v) = 1./r(v);
                        X = poly(r);
                        XX(g(i)+1:g(i+1)) = X(2:end);
                    end
            end
        end
    end

    function XX = bounceN(XX)
        for i = 1 : n
            switch t(i)
                case 1 % covariance
                    %%%% TODO: bounce negative covariance?
                    % warning('ssm:ssparam:NegCov', 'Negative covariance.\n');
                case 2 % ar>=3
                    for j = 1 : N
                        r   = roots([1 -XX(j, g(i)+1:g(i+1))]);
                        v   = abs(r) > 1;
                        if any(v)
                            r(v)    = 1./r(v);
                            X       = poly(r);
                            XX(j, g(i)+1:g(i+1)) = -X(2:end);
                        end
                    end
                case 3 % ma>=3
                    for j = 1 : N
                        r   = roots([1 XX(j, g(i)+1:g(i+1))]);
                        v   = abs(r) > 1;
                        if any(v)
                            r(v)    = 1./r(v);
                            X       = poly(r);
                            XX(j, g(i)+1:g(i+1)) = X(2:end);
                        end
                    end
            end
        end
    end
end

