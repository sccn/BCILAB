function B = subsref(A, S)

%@SSMODEL/SUBSREF Subscripted reference.
%   .name returns model name.
%   .Hinfo returns noise model information.
%   .info returns component information.
%   .p returns model data dimension.
%   .m returns model state dimension.
%   .r returns model state disturbance dimension.
%   .n returns model time duration.
%   .ncom returns number of model components.
%   .H returns model observation disturbance variance.
%   .Z returns model state to observation transformation.
%   .T returns model state transition.
%   .R returns model state disturbance to state transformation.
%   .Q returns model state disturbance variance.
%   .c returns model state transition constant.
%   .a1 returns model initial state vector.
%   .P1 returns model initial state variance.
%   .sta returns true if the model is stationary.
%   .linear returns true if both Z and T are linear.
%   .gauss returns true if both H and Q are Gaussian.
%   .w returns the number of parameters.
%   .paramname returns model parameter names.
%   .param returns model parameter values.
%   .psi returns transformed model parameter values.
%   .q returns number of diffuse elements.

% (c) 2006-2007 Jyh-Ying Peng ´^´¼·­
% $Revision 1.0.0 $ $Generated: 2007/09/04 $

switch S(1).type
    case '.'
        switch S(1).subs
            case 'name',        B = A.name;
            case 'Hinfo',       B = A.Hinfo;
            case 'info',        B = A.info;
            case 'p',           B = A.p;
            case 'm',           B = size(A.T, 1);
            case 'r',           B = size(A.R, 2);
            case 'n',           B = A.n;
            case 'ncom',        B = length(A.mcom) - 1;
            case 'H',           B = A.H;
            case 'Z',           B = A.Z;
            case 'T',           B = A.T;
            case 'R',           B = A.R;
            case 'Q',           B = A.Q;
            case 'c',           B = A.c;
            case 'a1',          B = A.a1;
            case 'P1',          B = A.P1;
            case 'sta',         B = issta(A);
            case 'linear',      B = islinear(A);
            case 'gauss',       B = isgauss(A);
            case 'w',           B = A.psi.w;
            case 'paramname',   B = A.psi.name;
            case 'param',       B = get(A.psi);
            case 'psi',         B = A.psi.value;
            case 'q',           B = nnz(A.P1.mat == Inf);
            otherwise,          error('ssm:ssmodel:subsref', 'Invalid subscripted reference.');
        end
    otherwise
        error('ssm:ssmodel:subsref', 'Invalid subscripted reference.');
end
if length(S) > 1, B = subsref(B, S(2:end)); end

