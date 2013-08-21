function [U,Out] = RecPF(m,n,aTV,aL1,picks,B,PsiT,Psi,opts,varargin,mfft2,imfft2,tvnesta,delta,bbb,AAA)
% [U,Out] = RecPF(m,n,aTV,aL1,picks,B,PsiT,Psi,opts,varargin)
%
% RecPF solves the TVL1-L2 model:
%
%   min aTV*TV(u) + aL1*|PsiT*U|_1 + 0.5|Fp*U - B|_2^2
%
% Inputs:
%
%  m, n     -- size of image
%  aTV, aL1 -- regularization parameters in the model
%  picks    -- sample positions in Fourier domain
%  B        -- measurment vector
%  PsiT     -- sparsifying basis, PsiT*U is the sparse representation of U
%  Psi      -- inverse of PsiT, Psi*x = U reconstructs the image
%  opts      --- contains parameters for algorithm
%              * opts.mit_inn: maxium inner iteration number for each beta {default 30}
%              * opts.mit_out: maxium outer iteration number {default 10} 
%              * opts.tol_inn: inner iteration error tolerance {1.e-3}
%                (this option is effective only when "opts.stc = 2", see below)
%              * opts.tol_rel_inn: tolerance of relative change in an inner iteration {5.e-2}
%              * opts.tol_rel_out: tolerance of relative change in an outer iteration {1.e-2}
%                (the above two options are effective only when "opts.stc = 1", see below)
%              * opts.beta0:   initial penalty parameter {2^5}
%              * opts.beta_max: final penalty parameter  {2^15}
%              * opts.beta_rate: increase rate of beta
%              * opts.idisp: 0 or nonzero, display inner iteration info or not
%              * opts.recordf: 0 or 1, keep (or not) history of function values, TV and fidelity;
%              * opts.U0: starting poing {default "a least squares solution"} 
%              * opts.wbarflag: 1 (on) or 0 (off), controls wait bar.
%              * opts.stc: 1 (stopping criterion based on Relative Change) 
%                          or 2 (stopping criterion based on optimality conditions) {default 1}
%              * opts.TVtype: 1 or 2, corresponding to anisotropic/isotropic TV {default 2}
%  varargin{1}  -- a m x n matrix contains local weights of TV. Specifically, the local weigthed 
%                  TV is discretized as:
%                     TV(u) = sum_i w_i ||D_i u||. 
%                  If all w_i == 1, then it is just the isotropic discritization of normal TV. 
%                  We require all w_i > 0.
%  varargin{2}  -- true image (used to compute relative errors when it is present)
%
% Outputs:
%     U   --- reconsctructed image
%     Out --- a structrue contains
%             * Out.iter: total iteration number
%             * Out.Inner: inner iteration numbers for each beta
%             * Out.ftrue: function values at each iteration 
%             * Out.Fvalvscpu: function values decrease with respect to CPU time 
%             * Out.TVhist: total variation changes with iteration 
%             * Out.FIDhist: fidelity history
%              {the above four subfields are present only when opts.recordf = 1}
%             * Out.ITvsERR: relative error vs. iteration number 
%             * Out.CPUvsERR: CPU time vs. relative error 
%              {the above two subfields are present only when the original image is present}

%
% Yin Zhang, 11-10-2007
% CAAM, Rice University, Copyright (2007)
%
% Junfeng Yang, last modified on Feb. 2, 2009
%

global Ux Uy FU PsiTU
global Wx Wy
global Z PsiZ
global beta
global remain

if aTV == 0 && aL1 == 0; error('No regularization' ); end
if aTV <  0 || aL1 <  0; error('Regularization < 0'); end
if ~exist('opts','var'); opts = []; end

[mit_inn,mit_out,tol_inn,tol_rel_inn,tol_rel_out,beta0,beta_max, ...
    beta_rate,TVtype,idisp,recordf,U0,wbarflag,stc] = setopts(opts);

[Nomin1,Denom1,Denom2] = getC(picks,B,aTV,m,n);

remain = setdiff(1:m*n,picks);
if isempty(U0); U = zeros(m,n); else U = U0; clear U0; end

if ~isempty(varargin)
    Weights = varargin{1};
    [mm,nn] = size(Weights);
    if mm~=m || nn~= n
        Weights = ones(m,n);
    elseif find(Weights < 0,1)
        error('negative weights are not allowed!');
    end
else
    Weights = ones(m,n);
end
if length(varargin) == 2
    I = varargin{2};
    nrmI = norm(I(:));  rer = norm(U(:)-I(:))/nrmI;
    RER = zeros(500,2); RER(1,:) = [0,rer];
    CPUvsRER = zeros(500,2); CPUvsRER(1,:) = [0,rer];
end

if aL1 > 0; PsiTU = PsiT(U); end

if aTV > 0
    Ux = [diff(U,1,2), U(:,1) - U(:,n)];
    Uy = [diff(U,1,1); U(1,:) - U(m,:)];
end

if recordf
    fval_true = funcval(U,aTV,aL1,B,picks);
    FVAL = zeros(500,2);
    FVAL(1,:) = [0,fval_true];
    Fvalvscpu = zeros(500,2);
    Fvalvscpu(1,:) = [0,fval_true];
    TVhist = zeros(500,1);
    FIDhist = zeros(500,1);
end

% initialization 
Inn = zeros(100,1);
beta = beta0;
totalIter = 0;
Outiter = 0;
stopc = 0;
t0 = cputime;

if wbarflag == 1
    wbar = waitbar(0, 'RecPF is running, please wait ...');
end

%% Main loop
while ~stopc
    Outiter = Outiter + 1;
    if stc == 1; Uo = U; end
    
    Denom = Denom1;
    if aTV > 0; Denom = Denom + (aTV*beta)*Denom2; end
    if aL1 > 0; Denom = Denom + aL1*beta; end

    stopci = 0; Inniter = 0; 
    while ~stopci

        % ================================
        %  Begin Alternating Minimization
        % ----------------
        %   W-subprolem
        % ----------------
        if aTV > 0;
            switch TVtype
                case 1;   % anisotropic TV
                    Wx = sign(Ux).* max(abs(Ux)-Weights./beta,0);
                    Wy = sign(Uy).* max(abs(Uy)-Weights./beta,0);
                case 2;   % isotropic TV
                    V = sqrt(Ux.^2 + Uy.^2);
                    S = max(V - Weights./beta, 0);
                    S = S ./ max(V,eps);
                    Wx = S.*Ux; Wy = S.*Uy;
                    clear V S
                otherwise; error('TVtype must be 1 or 2');
            end
        end

        % ----------------
        %   Z-subprolem
        % ----------------
        if aL1 > 0;
            Z = sign(PsiTU).*max(abs(PsiTU)-1/beta,0);
            PsiZ = Psi(Z);
        end

        % ----------------
        %   U-subprolem
        % ----------------
        if stc == 1; Ui = U; end
        rhs = 0;
        if aTV > 0
            rhs = [Wx(:,end) - Wx(:, 1), -diff(Wx,1,2)];
            rhs = rhs + [Wy(end,:) - Wy(1, :); -diff(Wy,1,1)];
            rhs = (aTV*beta)*rhs;
        end
        if aL1 > 0
            rhs = rhs + (aL1*beta)*PsiZ;
        end

        Nomin = Nomin1 + mfft2(rhs);

        FU = Nomin./Denom;
        U = real(imfft2(FU));
        Inniter = Inniter + 1;
        totalIter = totalIter + 1;
        %
        %  End Alternating Minimization
        % ================================

        % -----------------------------------------
        % check if inner stopping criterion is met
        %
        if stc == 1
            if aTV > 0
                Ux = [diff(U,1,2), U(:,1) - U(:,n)];
                Uy = [diff(U,1,1); U(1,:) - U(m,:)];
            end
            if aL1 > 0; PsiTU = PsiT(U); end
            chg_inn = norm(Ui(:) - U(:))/norm(U(:));
            stopci = ((chg_inn < tol_rel_inn) || (Inniter >= mit_inn));
        elseif stc == 2
            [res_wn,res_wz,wzr,res_Zn,res_Zz,Zzr,res_u] = checkopt(U,aTV,aL1,PsiT,picks,B,0,Weights,mfft2,mifft2);
            res = [res_wn,res_wz,res_Zn,res_Zz,res_u];
            stopci = (max(res) < tol_inn || (Inniter >= mit_inn));
        else
            error('please specify a stopping criterion.');
        end    
        
        if recordf && totalIter < 500
            [fval_true,tvu,fidu] = funcval(U,aTV,aL1,B,picks);
            TVhist(totalIter) = tvu;
            FIDhist(totalIter) = fidu;
            FVAL(totalIter+1,:) = [totalIter,fval_true];
            Fvalvscpu(totalIter+1,:) = [cputime - t0, fval_true];
        end
        if exist('I','var') && totalIter < 500
            rer = norm(U(:)-I(:))/nrmI;
            RER(totalIter+1,:) = [totalIter,rer];
            CPUvsRER(totalIter+1,:) = [cputime - t0, rer];
        end
        if idisp && stc == 2
            fprintf('Iter: %d, res_ wn: %4.1e, wz: %4.1e, Zn: %4.1e, Zz: %4.1e, u: %4.1e, wzr %2.0f%%, Zzr %2.0f%%\n', ...
                totalIter, res_wn,res_wz,res_Zn,res_Zz,res_u,100*wzr,100*Zzr);
        elseif idisp && stc == 1
            fprintf('Iter: %d, chg_inn %4.2f\n',totalIter,chg_inn);
        end

    end % inner
    Inn(Outiter) = Inniter;

    % ------------------------------------------
    % check if outer stopping criterion is met
    %
    beta = beta_rate*beta;
    if stc == 1
        chg_out = norm(U(:)-Uo(:))/norm(U(:));
        stopc = ((chg_out < tol_rel_out) || (Outiter >= mit_out));
    elseif stc == 2
        stopc = ((beta > beta_max) || (Outiter >= mit_out));
    else
        error('please specify a stopping criterion.');
    end
    
    %%%%%%% CHANGE HERE
    
    tv1 = calctv(m,n,U);
    R = fft2(U);
    R = R(picks);
    err1 =norm(bbb-AAA(U(:)));
    if (tv1 < 1.05*tvnesta)&(err1 < 1.05*delta)
        disp(sprintf('Out with modified stopping criterion - TV < 1.05 TV_nesta and ||b - A x||_2 < 1.05 ||b - A x_nesta||_2 '));
        break
    end
    
    if wbarflag == 1
        waitbar(log2(beta)/(log2(beta_max)+1), wbar)
    end
    
end % outer
Out.iter = totalIter;
Out.Inner = Inn(1:Outiter);
if wbarflag == 1
    close(wbar);
end
if recordf
    Out.ftrue = FVAL(1:totalIter+1,:);
    Out.Fvalvscpu = Fvalvscpu(1:totalIter+1,:);
    Out.TVhist = TVhist(1:totalIter);
    Out.FIDhist = FIDhist(1:totalIter);
end
if exist('I','var')
    Out.ITvsRER  = RER(1:totalIter+1,:);
    Out.CPUvsRER = CPUvsRER(1:totalIter+1,:);
end

%% ----------------- SUBFUNCTION ---------------------------
function [res_wn,res_wz,wzr,res_Zn,res_Zz,Zzr,res_u] = ...
    checkopt(U,aTV,aL1,PsiT,picks,f,flag,Weights,mfft2,imfft2)
%
% This function checks the optimality of
%                   (              beta                 )
% min   aTV * sum_i (||(W_i)||_2 + ---- ||D_iu - W_i||^2)
%                   (               2                   )
%
%                   (         beta             )
%     + aL1 *       ( |Z|_1 + ---- |Z-PsiT*u|^2 )
%                   (           2              )
%
%     + .5 |F_p*u - f|^2
%
% where W_i = (Wx_i;Wy_i);
%
% U     --- current point
% f     --- measurment vector
% flag  --- 0 (do not check res_u), otherwise check res_u;
%           Generally, res_u << 0, so setting flag = 0 is fine.
% Weights -- local weights of TV
%

% Junfeng Yang, Aug. 08, 2008

global Ux Uy FU PsiTU
global Wx Wy
global beta
global remain
global Z PsiZ

[m,n] = size(U);
mn = m*n;

% res_wn, res_wz, wzr
if aTV > 0
    Ux = [diff(U,1,2), U(:,1) - U(:,n)];
    Uy = [diff(U,1,1); U(1,:) - U(m,:)];

    Dx = Wx - Ux;
    Dy = Wy - Uy;

    V = sqrt(Wx.^2 + Wy.^2);
    Iz = find(V == 0);

    if isempty(Iz)
        In = 1:mn;
    else
        In = find(V ~= 0);
        V(Iz) = 1;
    end

    V = V*beta;

    RWx = Dx + Weights.*Wx./V;
    RWy = Dy + Weights.*Wy./V;
    S = sqrt(RWx.^2 + RWy.^2);

    wzr = length(Iz)/mn;
    if isempty(In)
        res_wn = 0;
    else
        res_wn = max(S(In));
    end

    if isempty(Iz);
        res_wz = -1;
    else
        res_wz = max(S(Iz) - Weights(Iz)/beta);
    end

else
    res_wn = 0;
    res_wz = 0;
    wzr = 0;
end


% res_Zn, res_Zz, Zzr
if aL1 > 0
    PsiTU = PsiT(U);
    absPsiTU = abs(PsiTU);

    ZIn = find(Z ~= 0);
    if isempty(ZIn)
        res_Zn = 0;
    else
        R3 = (1/beta)*sign(Z(ZIn)) + Z(ZIn) - PsiTU(ZIn);
        res_Zn = max(abs(R3));
    end

    ZIz = find(Z == 0);
    if isempty(ZIz)
        res_Zz = -1;
    else
        R4 = absPsiTU(ZIz) - 1/beta;
        res_Zz = max(R4);
    end

    Zzr = length(ZIz)/m/n;

else
    PsiTU = 0;
    res_Zn = 0;
    res_Zz = -1;
    Zzr = 0;
end

% res_u
res_u = 0;
if flag == 0
    return;
else
    % res_u
    FU(picks) = -f + FU(picks);
    FU(remain) = 0;
    Ru = imfft2(FU);
    Ru = real(Ru);

    if aTV > 0
        Dxx = [Dx(:,1) - Dx(:,n), diff(Dx,1,2)];
        Dyy = [Dy(1,:) - Dy(m,:); diff(Dy,1,1)];
        Ru = Ru + (beta*aTV)*(Dxx + Dyy);
    end

    if aL1 > 0
        Ru = Ru + (aL1*beta)*(U - PsiZ);
    end
    res_u = max(abs(Ru(:)));
end
%% ------------- SUBFUNCTION ---------------
function [Nomin1,Denom1,Denom2] = getC(picks,B,aTV,m,n)

% compute fixed quantities
Nomin1 = zeros(m,n);
Nomin1(picks) = B;
Denom1 = zeros(m,n);
Denom1(picks) = 1;
Denom2 = 0;
if aTV > 0
    Denom2 = abs(psf2otf([1,-1],[m,n])).^2 + abs(psf2otf([1;-1],[m,n])).^2;
end
%% ----------  SUBFUNCTION ---------------
function [mit_inn,mit_out,tol_inn,tol_rel_inn,tol_rel_out,beta0,beta_max, ...
    beta_rate,TVtype,idisp,recordf,U0,wbarflag,stc] = setopts(opts)

% define default option fields
mit_inn = 30;
mit_out = 10;
tol_inn  = 1.e-3; %--- 
tol_rel_inn = 5.e-3; %---
tol_rel_out = 1.e-2; %---
beta0 = 2^5;
beta_max = 2^15;
beta_rate = 2;
TVtype = 2;
idisp = 0;
recordf = 0;
U0 = [];
wbarflag = 0;
stc = 1;

% change to specified option fields if exist
if ~isempty(opts);
    if ~isa(opts,'struct'); error('L1pfi: opts not a struct'); end
    if isfield(opts,'mit_inn'); mit_inn = opts.mit_inn; end
    if isfield(opts,'mit_out'); mit_out = opts.mit_out; end
    if isfield(opts,'tol_inn'); tol_inn = opts.tol_inn; end
    if isfield(opts,'tol_rel_inn'); tol_rel_inn = opts.tol_rel_inn; end
    if isfield(opts,'tol_rel_out'); tol_rel_out = opts.tol_rel_out; end
    if isfield(opts,'beta0');   beta0 = opts.beta0; end
    if isfield(opts,'beta_max'); beta_max = opts.beta_max; end
    if isfield(opts,'beta_rate'); beta_rate = opts.beta_rate; end
    if isfield(opts,'TVtype'); TVtype = opts.TVtype; end
    if isfield(opts,'idisp'); idisp = opts.idisp; end
    if isfield(opts,'recordf'); recordf = opts.recordf; end
    if isfield(opts,'U0'); U0 = opts.U0; end
    if isfield(opts,'wbarflag'); wbarflag = opts.wbarflag; end
    if isfield(opts,'stc'); stc = opts.stc; end
end
