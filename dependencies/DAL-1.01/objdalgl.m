% objdalgl - objective function of DAL with grouped L1 regularization
%
% Copyright(c) 2009 Ryota Tomioka
% This software is distributed under the MIT license. See license.txt

function varargout=objdalgl(aa, info, prob, ww, uu, A, AT, B, lambda, eta)

nn=sum(info.blks);

vv   = ww/eta+AT(aa);
[vsth,ss] = prob.softth(vv, lambda, info);
nm=ss+lambda; %% only correct for active components
I = find(ss>0);

if nargout<=3
  [floss, gloss, hmin]=prob.floss.d(aa, prob.floss.args{:});
else
  [floss, gloss, hloss,hmin]=prob.floss.d(aa, prob.floss.args{:});
end

m = size(gloss,1);

fval = floss+0.5*eta*sum(ss.^2);

if ~isempty(uu)
  u1   = uu/eta+B'*aa;
  fval = fval + 0.5*eta*sum(u1.^2);
end

varargout{1}=fval;

if nargout<=2
  varargout{2}=info;
else
  gg  = gloss+eta*A(vsth);
  soc = sum((vsth-ww/eta).^2);
  if ~isempty(uu)
    gg  = gg+eta*B*u1;
    soc = soc+sum((B'*aa).^2);
  end

  if soc>0
    info.ginfo = norm(gg)/(sqrt(eta*hmin*soc));
  else
    info.ginfo = inf;
  end
  varargout{2} = gg;

  if nargout==3
    varargout{3} = info;
  else
    switch(info.solver)
      case 'cg'
       prec = hloss;
       varargout{3} = struct('blks',info.blks,'hloss',hloss,...
                             'I',I,'vv',vv,'nm',nm,'prec',prec,'lambda',lambda,'B',B);
     case 'nt'
      if length(I)>0
        AF = zeros(m, sum(info.blks(I)));
        ix0=0;
        for kk=1:length(I)
          jj=I(kk);
          bsz=info.blks(jj)
          J=sum(info.blks(1:jj-1))+(1:bsz);
          vn=vv(J)/nm(jj);
          ff=sqrt(1-lambda/nm(jj));

          Iout=ix0+(1:bsz);
          ix0=Iout(end);
          F1 = sparse(J,1:bsz,ff*ones(bsz,1),nn,bsz);
          F2 = sparse(J,ones(bsz,1),(1-ff)*vn,nn,1);
          AF(:,Iout)=A(F1)+A(F2)*vn';
        end
        H = hloss+ eta*AF*AF';
      end
      if ~isempty(uu)
        H = H+eta*B*B';
      end

      varargout{3} = H;
     case 'ntsv'
      H=hloss;
      for kk=1:length(I)
        jj=I(kk);
        bsz=info.blks(jj)
        J=sum(info.blks(1:jj-1))+(1:bsz);
        vn=vv(:,J)/nm(jj);
        ff=sqrt(1-lambda/nm(jj));

        F1 = sparse(J,1:bsz,ff*ones(bsz,1),nn,bsz);
        F2 = sparse(J,ones(bsz,1),(1-ff)*vn,nn,1);
        AF=A(F1)+A(F2)*vn';
        H = H+eta*AF*AF';
      end
      if ~isempty(uu)
        H = H+eta*B*B';
      end
      varargout{3} = H;
    end % end switch(info.solver)
    varargout{4} = info;
  end
end

