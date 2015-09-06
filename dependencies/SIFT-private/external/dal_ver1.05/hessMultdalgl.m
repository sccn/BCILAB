

% hessMultdalgl - function that computes H*x for DAL with grouped
%                 L1 regularization
%
% Copyright(c) 2009 Ryota Tomioka
% This software is distributed under the MIT license. See license.txt

function yy = hessMultdalgl(xx, A, eta, Hinfo)

yy = Hinfo.hloss*xx;
if ~isempty(Hinfo.precomp)
    % general case
    for p=Hinfo.precomp
        % project & reshape xx's for all spans
        xk = reshape(p.AJ'*xx,[],numel(p.jj));
        % dot product between vn's and xk's, scale by ff's and multiply by vn's (also add xk*(1-ff))
        tmp = bsxfun(@times,p.vn,sum(p.vn.*xk).*p.ff) + bsxfun(@times,xk,p.omff);
        % block-diagonalize tmp
        bd = sparse(p.idxu,p.idxv,tmp,numel(xk),numel(p.jj));
        % map through AJ again, sum and multiply by eta(1), add results to yy
        yy = yy + sum(p.AJ*bd,2)*eta(1);
    end
else
    % very sparse case
    blks =Hinfo.blks;
    hloss=Hinfo.hloss;
    I    =Hinfo.I;
    vv   =Hinfo.vv;
    nm   =Hinfo.nm;
    lambda=Hinfo.lambda;
    for kk=1:length(I)
        jj=I(kk);
        J=Hinfo.blkival{jj};
        vn=vv(J)/nm(jj);        
        ff=lambda/nm(jj);
        AJ=A.slice(J);
        xk=AJ'*xx;
        yy = yy + eta(1)*(AJ*((1-ff)*xk + ff*(vn'*xk)*vn));
    end
end

B=Hinfo.B;
if ~isempty(B)
  yy = yy + eta(2)*(B*(B'*xx));
end