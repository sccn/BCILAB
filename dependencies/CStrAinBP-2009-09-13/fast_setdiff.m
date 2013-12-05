function [A,A_in_B] = fast_setdiff(A,B)
% A fast version of setdiff for cell arrays of strings.

try
    if isempty(A)
        A = {}; 
    else
        A_in_B = CStrAinBP(A,B);
        A(A_in_B) = [];
    end
catch %#ok<CTCH>
    disp_once('Using the slower MATLAB fallback for setdiff(). Consider setting up a compiler to get much improved BCILAB performance (in some areas).');
    if nargout <= 1
        A = setdiff(A,B);
    else
        A_in_B = 1:length(A);
        [A,A_not_in_B] = setdiff(A,B);
        A_in_B(A_not_in_B) = [];
    end
end
