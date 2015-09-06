function [outerproduct] = outer(tensor1, tensor2, squeezedimensions)
%C = OUTER(A, B, squeezedimensions) Computes the outer product
%C(i[1],...,i[m],j[1],...,j[n]) = A(i[1],...,i[m]) * B(j[1],...j[n]).
%Discards superfluous singleton dimensions if squeezedimensions ~= 0.
%
%Note: thanks to Emese Toth for suggesting the algorithm used here.
%
%Wynton Moore, January 2006


%store initial dimensions
dim1=size(tensor1);dim2=size(tensor2);


%evaluate
outerproduct=reshape(reshape(tensor1, [], 1)*reshape(tensor2, 1, []), [dim1 dim2]);


%discard superfluous singleton dimensions
if squeezedimensions
    outerproduct=squeeze(outerproduct);
end
