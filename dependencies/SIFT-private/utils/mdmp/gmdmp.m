function [gmdmproduct] = gmdmp(tensor1, d1, tensor2, d2)
%GMDMP General Multi Dimensional Matrix Product.
%C = GMDMP(A, d1, B, d2) Computes the product
%C(i[1],...,i[d1-1],i[d1+1],...,i[m],j[1],...,j[d2-1],j[d2+1],...,j[n]) =
%     A(i[1],...,i[d1-1],k,i[d1+1],...,i[m]) * B(j[1],...,j[d2-1],k,j[d2+1],...,j[n])
%(Sum on k).
%
%C = GMDMP(A, d1, B, d2) takes the outer product of A and B, then traces
%along the diagonal formed by dimensions d1 of A and d2 of B. For example,
%C = GMDMP(A, ndims(A), B, 1) is just the natural extension of 2D matrix
%multiplication, and for A and B both 2D, coincides with C = A * B.
%
%Note: it matters not if the lengths of dimensions d1 of A and d2 of B do
%not agree (see header of DIAGSUM.M).
%
%Wynton Moore, January 2006


%evaluate
gmdmproduct=diagsum(outer(tensor1, tensor2, 0), d1, ndims(tensor1)+d2);
