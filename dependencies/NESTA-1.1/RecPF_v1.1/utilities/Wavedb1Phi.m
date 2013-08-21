function Y = Wavedb1Phi(X,trans)
persistent s;

level = 2;

if ~trans;
    % Phi' * X
    [Y,s] = wavedec2(X,level,'db1');
else
    Y = waverec2(X(:),s,'db1');
end
    
 
