% objdalgl - objective function of DAL with grouped L1 regularization
%
% Copyright(c) 2009 Ryota Tomioka
% This software is distributed under the MIT license. See license.txt

function [fval,out2,out3,out4] = objdalgl(aa, info, prob, ww, uu, A, B, lambda, eta)

% A' * aa
if isempty(info.ATaa)
    info.ATaa = A.Ttimes(aa);  end

% soft-thresholding
vv = ww/eta(1)  + info.ATaa;
[vsth,ss] = prob.softth(vv, lambda, info);

% new weights and spectrum
info.wnew = vsth*eta(1);
info.spec = ss;

switch nargout
    case 3
        % implicit-Hessian methods (qn)        
        [floss, gloss, hmin] = prob.floss.d(aa, prob.floss.args{:});
        gg  = gloss + eta(1)*A.times(vsth);
        soc = sum((vsth - ww/eta(1)).^2);
        
        % compute objective-function value
        if isempty(uu)
            % no unregularized components
            fval = floss + 0.5*eta(1)*sum(ss.^2);
        else
            % have uregularized components
            u1   = uu/eta(2)+B'*aa;
            gg  = gg + eta(2)*B*u1;
            soc = soc + sum((B'*aa).^2);
            fval = floss + 0.5*eta(1)*sum(ss.^2) + 0.5*eta(2)*sum(u1.^2);
        end
        
        info.ginfo = norm(gg) / (sqrt(min(eta) * hmin * soc));
        out2 = gg;
        out3 = info;
        
    case 4
        % explicit-Hessian methods (cg,nt,ntsv)
        [floss, gloss, hloss, hmin] = prob.floss.d(aa, prob.floss.args{:});
        gg  = gloss + eta(1)*A.times(vsth);
        soc = sum((vsth - ww/eta(1)).^2);

        nm = ss+lambda;
        nn = sum(info.blks);
        m = size(gloss,1);
        
        % compute objective-function value
        if isempty(uu)
            % no unregularized components
            fval = floss + 0.5*eta(1)*sum(ss.^2);
        else
            % have uregularized components
            u1   = uu/eta(2)+B'*aa;
            gg  = gg + eta(2)*B*u1;
            soc = soc + sum((B'*aa).^2);
            fval = floss + 0.5*eta(1)*sum(ss.^2) + 0.5*eta(2)*sum(u1.^2);
        end
        
        out2 = gg;
        info.ginfo = norm(gg) / (sqrt(min(eta) * hmin * soc));
        
        % non-zero indices
        I = find(ss>0);
        
        switch(info.solver)
            case 'cg'
                prec = hloss;
                if ~isempty(uu)
                    prec = prec+spdiag(eta(2)*sum(B.^2,2)); end
                % pre-compute a few things for the Hessian product
                k=1;
                if length(I) > 3
                    % for each set of non-zero indices in same-sized blocks...
                    for g=1:length(info.blksizes)
                        % get the indices
                        jj = I(info.blks(I)==info.blksizes(g))';
                        if ~isempty(jj)
                            % concatenate their block index spans
                            J = [info.blkival{jj}];
                            % read the non-zeros from inv(nm)
                            inm = 1./ nm(jj)';
                            % read out the joint block span from vv, reshape into [blocklen x #indices], normalize
                            vn = bsxfun(@times,reshape(vv(J),[],numel(jj)),inm);
                            idxu = 1:numel(vn); idxv = floor(1:1/size(vn,1):numel(jj)+1-1/size(vn,1));
                            % get ff's
                            ff=lambda*inm;
                            omff = 1-ff;
                            % read out the spans from A
                            AJ = A.slice(J);
                            % AJ, jj, vn, ff, omff, idxu, idxv
                            precomp(k) = struct('AJ',AJ,'jj',jj,'vn',vn,'ff',ff,'omff',omff,'idxu',idxu,'idxv',idxv);
                            k=k+1;
                        end
                    end
                else
                    precomp = [];
                end
                % dump a whole bunch of stuff in the output...
                out3 = struct('precomp',precomp,'blks',info.blks,'blkvec',{info.blkvec},'blkgrp',{info.blkgrp},'blkival',{info.blkival},'blksizes',info.blksizes,'hloss',hloss,...
                    'I',I,'vv',vv,'nm',nm,'prec',prec,'lambda',lambda,'B',B);
            case 'nt'
                % warning: horribly unoptimized
                H = hloss;
                if length(I)>0
                    AF = zeros(m, sum(info.blks(I)));
                    ix0  = 0;
                    cumblks = [0,cumsum(info.blks)];
                    for kk=1:length(I)
                        jj = I(kk);
                        bsz = info.blks(jj);
                        J = cumblks(jj)+(1:bsz);
                        vn = vv(J)/nm(jj);
                        ff = sqrt(1-lambda/nm(jj));
                        
                        Iout = ix0+(1:bsz);
                        ix0 = Iout(end);
                        F1 = sparse(J,1:bsz,ff*ones(bsz,1),nn,bsz);
                        F2 = sparse(J,ones(bsz,1),(1-ff)*vn,nn,1);
                        AF(:,Iout) = A.times(F1)+A.times(F2)*vn';
                    end
                    H = H + eta(1)*AF*AF';
                end
                if ~isempty(uu)
                    H = H+eta(2)*B*B'; end
                out3 = H;
            case 'ntsv'
                % warning: horribly unoptimized
                H = hloss;
                for kk=1:length(I)
                    jj = I(kk);
                    bsz = info.blks(jj)
                    J = sum(info.blks(1:jj-1))+(1:bsz);
                    vn = vv(:,J)/nm(jj);
                    ff = sqrt(1-lambda/nm(jj));
                    F1 = sparse(J,1:bsz,ff*ones(bsz,1),nn,bsz);
                    F2 = sparse(J,ones(bsz,1),(1-ff)*vn,nn,1);
                    AF = A.times(F1)+A.times(F2)*vn';
                    H = H+eta(1)*AF*AF';
                end
                if ~isempty(uu)
                    H = H+eta(2)*B*B'; end
                out3 = H;
        end
        out4 = info;
        
    case {1,2}
        % other stuff
        fval = prob.floss.d(aa, prob.floss.args{:}) + 0.5*eta(1)*sum(ss.^2);
        % add unregularized components
        if ~isempty(uu)
            u1   = uu/eta(2)+B'*aa;
            fval = fval + 0.5*eta(2)*sum(u1.^2);
        end
        out2 = info;
        return;
end
