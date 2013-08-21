function Y = dctPhi(X,trans)
if ~trans;
    Y = dct2(X);
else
    Y = idct2(X);
end
