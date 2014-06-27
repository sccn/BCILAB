classdef ParadigmPAL < ParadigmBaseSimplified
    % Experimental implementation of the Pattern Alignment Learning (PAL) method.
    %
    % This paradigm is not yet completed, please move on! :-)
    %
    % Name:
    %   Wave Alignment, work in progress
    %
    %                            Christian Kothe, Swartz Center for Computational Neuroscience, UCSD
    %                            2011-08-28
    
    methods
        
        function defaults = preprocessing_defaults(self)
            % define the default pre-processing parameters of this paradigm
            defaults = {'SpectralSelection',[0.1 35],'EpochExtraction',[-0.5 1.5], 'Resampling',100};
        end
        
        function model = calibrate_simple(self,varargin)
            % Override this function to implement your calibration code
            % Model = calibrate_simple(Signal,Arguments...)
            %
            % In:
            %   Signal : a single continuous EEGLAB data set (annotated with target markers; see
            %            set_targetmarkers for more info)
            %
            %   Arguments : further optional user arguments
            %
            % Out:
            %   Model : a model struct with a mandatory field .filter_graph and arbitrary other content
            %           * The .filter_graph filed is a 1x1 cell array that contains the desciption of
            %             filter steps that is to be applied to the data, and is usually passed as either
            %             the .tracking.online_expression field of the processed Signal or the processed
            %             Signal itself. If no signal processing is performed by this paradigm, the raw
            %             signal may be passed.
            %
            %           * May have optional fields .prediction_function, .prediction_window and
            %             .prediction_channels - though these are generally auto-deduced
            %
            % Notes:
            
            args = arg_define(varargin, ...
                arg_norep('signal'), ...
                arg_sub({'flt','SignalProcessing'}, self.preprocessing_defaults(), @flt_pipeline, 'Signal processing stages. These parameters control filter stages that run on the signal level; they can be enabled, disabled and configured for the given paradigm. The prediction operates on the outputs of this stage.'), ...
                arg({'lambdas','Lambdas'},[2.^(8:-0.125:-5)],[0 Inf],'Regularization parameter. Larger means stronger regularization, 1.0 is the Bayesian ad hoc choice.'), ...
                arg({'duration','WaveDuration'},0.8,[0 Inf],'Duration of Waveform. This is the length of the waveform that will be searched in the epoch.'), ...
                arg({'stepping','Stepping'},3,uint32([1 100]),'Relative jitter stepping. This is the precision (in samples) at which relative jitters will be considered.'), ...
                arg({'admm_iters','ADMMIterations'},30,uint32([1 1000]),'Max Iterations. Number of iterations for the overall ADMM solver.'), ...
                arg({'inner_iters','InnerIterations'},10,uint32([1 1000]),'Max Inner Iterations. Number of iterations for the inner ADMM solver that solves the sparsity vs trace norm.'), ...
                arg({'cg_iters','CgIterations'},30,uint32([1 1000]),'Max CG Iterations. Number of iterations for the inner conjugate-gradient solver.'), ...
                arg({'lanczos_iters','LanczosIterations'},10,uint32([1 1000]),'Max Lanczos Iterations. Number of iterations for the inner Lanczos solver.'), ...
                arg({'abs_tol','AbsoluteTolerance'},1e-6,[],'Absolute tolerance.'), ...
                arg({'rel_tol','RelativeTolerance'},1e-3,[],'Relative tolerance.'), ...
                arg({'alpha','OverRelaxation'},1,[1 1.8],'Over-relaxation parameter. Some references suggest that values between 1.5 and 1.8 can improve convergence.'), ...
                arg({'rho','AugmentedLagrangian'},1.0,[0 Inf],'Augmented Lagrangian parameter.'), ...
                arg({'rho_update','RhoUpdate'},true,[],'Update Rho. Whether to update rho dynamically according to 3.4.1 in [1]'), ...
                arg({'rho_cutoff','RhoUpdateThreshold'},10.0,[0 Inf],'Rho update threshold.'), ...
                arg({'rho_incr','RhoUpdateIncr'},2.0,[1 Inf],'Rho update increment factor.'), ...
                arg({'rho_decr','RhoUpdateDecr'},2.0,[1 Inf],'Rho update decrement factor.'), ...                
                arg({'kappa','SparsityTraceTradeoff'},0.75,[0 1],'Sparsity (=1) vs. Trace norm (=0) tradeoff. Between 0 and 1, similar to the elastic net criterion.'), ...
                arg({'patterns','NumPatterns'},20,uint32([1 1000]),'Number of weight patterns. This is the maximum number of (possibly distinct) weight patterns to recover.'), ...
                arg({'verbose','Verbose'},'really',[],'Verbose output.'), ...
                arg({'scaling','Scaling'}, 'std', {'none','center','std','minmax','whiten'}, 'Pre-scaling of the data. For the regulariation to work best, the features should either be naturally scaled well, or be artificially scaled.'), ...
                arg({'fast','Fast'},true,[],'Fast mode.'), ...
                arg({'arg_dialogsel','ConfigLayout'},self.dialog_layout_defaults(),[],'Parameters displayed in the config dialog. Cell array of parameter names to display (dot-notation allowed); blanks are translated into empty rows in the dialog. Referring to a structure argument lists all parameters of that struture, except if it is a switchable structure - in this case, a pulldown menu with switch options is displayed.','type','cellstr','shape','row'));
            
            % first pre-process the data (symbolically)
            signal = flt_pipeline('signal',args.signal, args.flt); %#ok<*NODEF>
            
            % evaluate this in an optimized fashion
            signal = exp_eval_optimized(signal);
            
            % get trial labels
            labels = set_gettarget(signal);
            
            % build meta-data
            model.args = rmfield(args,'signal');
            model.classes = unique(labels,'rows');
            model.tracking.filter_graph = signal;
            
            if length(model.classes) > 2
                error('Currently, only two classes are supported here.'); end
            
            % extract time-shifted windows from each trial...
            duration = round(args.duration * signal.srate);
            offsets = 1 : args.stepping : (signal.pnts - duration);
            trials = cell(1,length(offsets));

            % get the entire epoch and figure out an appropriate rescaling
            alldata = reshape(permute(signal.data,[2 1 3]),[],signal.trials);
            sc_info = hlp_findscaling(alldata,args.scaling);
            
            % for each time shift...
            for o = offsets
                % determine the window to cut out of the data
                cutwindow = o + (1:duration);
                % extract from all trials
                data = signal.data(:,cutwindow,:);
                % vectorize into #trials x #features, and rescale
                data = permute(data,[2 1 3]);
                trials{o} = hlp_applyscaling(reshape(data,[],signal.trials),sc_info)';
            end
            
            T = signal.trials;
            fprintf('#trials=%.0f; #shifts=%.0f; #features=%.0f; #channels=%.0f; #windows=%.0f\n',T,length(offsets),signal.nbchan*length(cutwindow),signal.nbchan,length(cutwindow));
            
            % vertically concatenate all shifts (into #shifts*#trials x #features)
            A = vertcat(trials{:});

            % remap the labels to -1 / +1
            labels(labels == model.classes(1)) = -1;
            labels(labels == model.classes(2)) = +1;
            
            % and generate shift-replicated target matrix b
            b = repmat(labels,1,length(offsets));
            b = b(:);
                
            C0 = -[b A];   % (bias is first elem, followed by rest...)
            ccp = C0';
            
            % weight vectors
            [m,n] = size(A);
            Wz = zeros(n+1,m);          % consensus variable
            Wx = zeros(n+1,m);          % data weights variable, FxN; n+1 accounts for the added bias (this is actually just the initial conditions vector...)
            Wu = zeros(n+1,m);          % scaled dual parameter u for Wx

            alpha = args.alpha;         % over-relaxation parameter (if 1, OR turned off)
            rho = args.rho;             % augmented Lagrangian parameter; probably determines sort of the tightness of the constraints...
            kappa = args.kappa;         % blend factor between trialwise sparsity and trialwise agreement
            K = args.patterns;          % the maximum number of SV's to recover
            
            rescale = 1;                  % scale of the left-hand side term (could be m/length(offsets)
            for mu = args.lambdas       % regularization parameter (should be decreasing)
                fprintf('Now using lambda = %.3f\n',mu);
                if args.verbose
                    fprintf('%3s\t%10s\t%10s\t%10s\t%10s\t%10s\n', 'iter', ...
                        'r norm', 'eps pri', 's norm', 'eps dual', 'objective');
                end
                
                for k = 1:args.admm_iters
                    Wzold = Wz;
                    fprintf('Outer step #%.0f\n',k);
                    
                    %% x-update
                    Wx = ParadigmPAL.minimize_batch(Wx, @ParadigmPAL.x_objective_batch, args.cg_iters, Wz, Wu, ccp, rho);
                    Wx_hat = alpha*Wx + (1-alpha)*Wzold;  % (optional over-relaxation)
                    
                    %% z-update
                    Wz = Wx_hat + Wu;
                    [Wz,fuuu,U,V] = ParadigmPAL.prox_sum(Wz,rescale,mu,kappa,rho,K,T,m,args);

                    %% dual parameter update
                    fprintf('  updating u''s\n');
                    Wu = Wu + (Wx_hat - Wz);
                    
                    %% diagnostics, reporting, termination checks
                    fprintf('  diagnostics...\n');
                    history.objval(k)  = 0; %ParadigmPAL.full_objective(ccp, mu, Wx, Wz);
                    history.r_norm(k)  = norm(Wx(:) - Wz(:));
                    history.s_norm(k)  = norm(rho*(Wz(:) - Wzold(:)));
                    history.eps_pri(k) = args.abs_tol + args.rel_tol*max(norm(Wx(:)), norm(Wz(:)));
                    history.eps_dual(k)= args.abs_tol + args.rel_tol*norm(rho*Wu(:));
                    
                    if args.verbose
                        fprintf('%3d\t%10.4f\t%10.4f\t%10.4f\t%10.4f\t%10.2f\n', k, ...
                            history.r_norm(k), history.eps_pri(k), ...
                            history.s_norm(k), history.eps_dual(k), history.objval(k));
                        if strcmp(args.verbose,'really') && nnz(fuuu)
                            f=figure; set(f,'Position',get(f,'Position')-[0 500 0 0]);
                            imagesc(Wz); ylabel('features'); xlabel('shifted trials'); title(sprintf('ADMM iter #%.0f, lambda=%.3f',k,mu)); 
                            figure;
                            % weight patterns                            
                            N = nnz(fuuu);
                            cols = ceil(sqrt(N));
                            rows = ceil(N/cols);
                            for p=1:N
                                subplot(cols,rows,p);
                                plot((0:round(args.duration * signal.srate)-1)/signal.srate,reshape(U(2:end,p),length(cutwindow),signal.nbchan)); title(sprintf('weight array #%.0f',p)); xlabel('time (s)'); ylabel('weight');
                            end
                            % topoplots...
                            f=figure; set(f,'Position',get(f,'Position')-[550 0 0 0]);
                            N = nnz(fuuu);
                            cols = ceil(sqrt(N));
                            rows = ceil(N/cols);
                            for p=1:N
                                subplot(cols,rows,p);
                                mat = reshape(U(2:end,p),length(cutwindow),signal.nbchan);
                                mat_inner = mat(round(end*1/4):round(end*3/4),:);
                                [dummy,peakidx] = max(mean(abs(mat_inner),2)); % #ok<NASGU>
                                topoplot(mat_inner(peakidx,:),signal.chanlocs);
                            end
                            drawnow; 
                        end
                    end
                    
                    if history.r_norm(k) < history.eps_pri(k) && history.s_norm(k) < history.eps_dual(k)
                        break; end
                    
                    %% update of rho
                    if args.rho_update
                        if history.r_norm(k) > args.rho_cutoff * history.s_norm(k)
                            disp('  incrementing rho.');
                            rho = rho .* args.rho_incr;
                            Wu = Wu ./ args.rho_incr;
                        elseif history.s_norm(k) > args.rho_cutoff * history.r_norm(k)
                            disp('  decrementing rho.');
                            rho = rho ./ args.rho_incr;
                            Wu = Wu .* args.rho_incr;
                        end
                    end
                end
            end
            1       
            %% plotting...
            
%             % full weight matrix
%             figure;imagesc(Wz); title('full weight matrix'); xlabel('features'); ylabel('shifted trials');
%             
%             % shift matrix (?)
%             figure; imagesc(reshape(V(:,1),length(offsets),[])); title('shift matrix'); xlabel('trials'); ylabel('shifts');
        end       
        
        function outputs = predict_simple(self,signal,model)
            % Override this function to implement your prediction code
            % Outputs = predict_simple(Sginal,Model)
            %
            % In:
            %   Signal : a signal pre-processed according to the model's filter graph
            %
            %   Model : a predictive model as created by your calibrate_simple() function
            %
            % Out:
            %   Outputs : a prediction/estimate for the most recent time point in the data (or one for
            %             every epoch if the signal is epoched); see ml_predict for the allowed formats
            
            error('This paradigm implements no predict() function.');
        end
        
        
        function visualize(self,model)
            % Optionally override this function to implement your visualization code
            % visualize(Model)
            %
            % In:
            %   Model: a model as created by your calibrate() function;
            %          a plot or GUI will be produced to inspect the model
            %
            
            error('This paradigm implements no visualize() function.');
        end
        
        
        
        % --- implementation ----
        
        function model = calibrate(self,varargin)
            % Implementation of the calibrate() function (see ParadigmBase)
            % Model = calibrate(Collection,GoalIdentifier,Arguments...)
            %
            % This is essentially a wrapper around calibrate_simple() which peels off any additional
            % data structure (for multiple collections and stream bundles) before passing it on.
            %
            % In:
            %   Collection : a collection (cell array) of stream bundles
            %
            %   GoalIdentifier: the goal-identifier
            %
            %   Arguments... : further optional user arguments
            %
            % Out:
            %   Model : The returned model struct
            
            % get the collection argument and extract the signal argument
            collection = arg_extract(varargin,{'collection','Collection'},[],{});
            if ~isempty(collection)
                if length(collection) > 1
                    error('This paradigm does not support collections of more than one recording.'); end
                if length(collection{1}.streams) > 1
                    error('This paradigm does not support more than one stream in parallel.'); end
                signal = collection{1}.streams{1};
            else
                signal = {};
            end
            
            % call calibrate_simple with the signal passed in as additional argument
            model = self.calibrate_simple('signal',signal,varargin{:});
        end
        
        
        function outputs = predict(self,bundle,model)
            % Override this function to implement your prediction code
            % Outputs = predict(Bundle,Model)
            %
            % This is a wrapper around the predict_simple() function which peels off the stream bundle
            % representation and just passes on the contained signal.
            %
            % In:
            %   Bundle : stream bundle preprocessed by the model's filter graph
            %
            %   Model : a predictive model
            %
            % Out:
            %   Outputs : the outputs of calibrate_simple()
            
            if length(bundle.streams) > 1
                error('This paradigm does not support more than one stream in parallel.'); end
            outputs = self.predict_simple(bundle.streams{1},model);
        end
        
    end
    
    
    methods(Static)
        function obj = full_objective(ccp, mu, Wx, Wz)
            % this is the overall objective function of alignment learning
            obj = sum(log(1 + exp(sum(ccp.*Wx)))) + mu*(norm_nuc(Wz) + sum(norms(Wz))) * size(ccp,2);            
        end
        
        function [val,grad] = x_objective(x, C, z, u, rho)
            % this is the objective function for the decoupled x criterion (basically l2 logreg)
            ecx = exp(C*x);
            val = sum(log(1 + ecx)) + (rho/2)*norm(x - z + u).^2;
            grad = C'*(ecx./(1 + ecx)) + rho*(x - z + u);
        end
        
        function [vals,grads] = x_objective_batch(Wx,Wz,Wu,ccp,rho)
            % this is the batch version of the above (for all samples & variables in parallel)
            ecx = exp(sum(ccp.*Wx));
            d = Wx - Wz + Wu;
            vals = log(1 + ecx) + (rho/2)*sum(d.*d);
            grads = bsxfun(@times,ccp,(ecx./(1 + ecx))) + rho*d;
        end
        
        % l1/l2 shrinkage operator
        function Wz = shrinkage(Wa, kappa)
            Wz = pos(1 - kappa/norm(Wa))*Wa;
        end

        % l1 shrinkage operator
        function z = shrinkage_l1(a, kappa)
            z = max(0, a-kappa) - max(0, -a-kappa);
        end

        % trace norm proximal operator
        function [W,fuuu,U,V] = prox_tr(W,rescale,mu,kappa,rho,K,args)
            fprintf('  trace norming r\n');
            [U,S,V] = pca(W,K,args.lanczos_iters);
            fuuu = ParadigmPAL.shrinkage_l1(diag(S),(rescale*mu*(1-kappa))/rho); nn = length(fuuu); % we shrink the SV spectrum using prox operator
            Ws = sparse(1:nn,1:nn,fuuu,size(S,1),size(S,2));
            % clamp minority of V's members at 0
            smv = sign(median(V));
            V(:,smv==-1) = min(0,V(:,smv==-1));
            V(:,smv==1) = max(0,V(:,smv==1));
            W = U*(Ws*V');
            
            if strcmp(args.verbose,'ultra')
                f=figure; set(f,'Position',get(f,'Position')-[0 500 0 0]);
                imagesc(W); ylabel('features'); xlabel('shifted trials'); title(sprintf('lambda=%.3f',mu));
                figure;
                % weight patterns
                N = nnz(spec);
                cols = ceil(sqrt(N));
                rows = ceil(N/cols);
                signal=args.signal;
                for p=1:N
                    subplot(cols,rows,p);
                    plot((0:(numel(U(2:end,p))/signal.nbchan)-1)/signal.srate,reshape(U(2:end,p),[],signal.nbchan)); title(sprintf('weight array #%.0f',p)); xlabel('time (s)'); ylabel('weight');
                end
                drawnow;
            end
            
        end
        
        % group-sparse proximal operator
        function W = prox_gl(W,rescale,mu,kappa,rho,T,m)
            fprintf('  sparsifying z\n');
            for i = 1:m
                W(2:end,i) = ParadigmPAL.shrinkage(W(2:end,i), (rescale*mu*kappa/10)/rho); end
        end
        
        % summed proximal operator (using some splitting scheme, see Bauscke, 2008)
        function [x,fuuu,U,V] = prox_sum(r,rescale,mu,kappa,rho,K,T,m,args)
            y = r;
            p = zeros(size(y));
            q = zeros(size(y));
            for n=1:args.inner_iters
                [x,fuuu,U,V] = ParadigmPAL.prox_tr(y+q,rescale,mu,kappa,rho,K,args);
                q = y+q-x;
                y = ParadigmPAL.prox_gl(x+p,rescale,mu,kappa,rho,T,m);
                p = x+p-y;
                if norm(x(:) - y(:)) < args.rel_tol*max(norm(x(:)),norm(y(:)))
                    break; end
            end            
        end
  
        function X = minimize_batch(X, f, len, Wz,Wu,ccp,rho)
            % Copyright (C) 2001 - 2010 by Carl Edward Rasmussen, 2010-01-03
            
            INT = 0.1;    % don't reevaluate within 0.1 of the limit of the current bracket         % const
            EXT = 3.0;                  % extrapolate maximum 3 times the current step-size         % const
            MAX = 20;                         % max 20 function evaluations per line search         % const
            RATIO = 10;                                       % maximum allowed slope ratio         % const
            SIG = 0.1; RHO = SIG/2; % SIG and RHO are the constants controlling the Wolfe-          % const
            if max(size(len)) == 2, red=len(2); len=len(1); else red=1; end                         % scalar

            
            [f0 df0] = feval(f,X,Wz,Wu,ccp,rho);          % get function value and gradient
            N = length(f0);
                        
            [A,B,d1,d2,d3,d4,f1,f2,f3,f4,x1,x2,x4] = deal(zeros(1,N)); % initialize a few variables...
            df3 = zeros(size(X));
            i = 0;                                            % zero the run length counter         % scalar
            ls_failed = false(1,N);                             % no previous line search has failed% boolean vector
            i = i + (len<0);                                            % count epochs?!            % scalar
            s = -df0;                                 % initial search direction (steepest)         % column MATRIX
            d0 = -sum(s.*s);                                                    % and slope         % row vector
            x3 = red./(1-d0);                                  % initial step is red/(|s|+1)        % row vector
            
            keep_going = true(1,N);                                                                 % mask for the outermost loop updates...            
            while any(keep_going) && i < abs(len)                                      % while not finished            % scalar
                i = i + (len>0);                                      % count iterations?!          % scalar
                fprintf('  Cg step #%.0f\n',i);
                
                X0 = X;                                                                             % column MATRIX
                F0 = f0;                                                                            % row vector
                dF0 = df0;                   % make a copy of current values                        % column MATRIX
                if len>0, M = MAX; else M = min(MAX, -len-i); end                                   % scalar
                
                om = keep_going;                                                                    % mask for outer-loop updates
                while 1                             % keep extrapolating as long as necessary
                    x2(om) = 0;                                                                     % row vector
                    f2(om) = f0(om);                                                                % row vector
                    d2(om) = d0(om);                                                                % row vector
                    f3(om) = f0(om);                                                                % row vector
                    df3(:,om) = df0(:,om);                                                          % column matrix
                    ok(om) = false;                                                                 % boolean vector
                    while ~all(ok(om)) && M > 0                                                     % mask check
                        M = M - 1; i = i + (len<0);                         % count epochs?!        % scalar
                        nokom = ~ok&om;
                        [f3(nokom) df3(:,nokom)] = feval(f, X(:,nokom) + bsxfun(@times,x3(nokom),s(:,nokom)), Wz(:,nokom), Wu(:,nokom), ccp(:,nokom), rho); % masked array op
                        ok(om) = ok(om) | ~(isnan(f3(om)) | isinf(f3(om)) | any(isnan(df3(om))+isinf(df3(om))));            % mask update
                        x3(~ok&om) = (x2(~ok&om)+x3(~ok&om))/2;         % bisect and try again      % maked update
                    end
                    mask = f3<F0;                                       % keep best values          % mask generation
                    X0(:,mask&om) = X(:,mask&om) + bsxfun(@times,x3(mask&om),s(:,mask&om));         % masked updates
                    F0(mask&om) = f3(mask&om);
                    dF0(:,mask&om) = df3(:,mask&om);
                    
                    d3(om) = sum(df3(:,om).*s(:,om));                                    % new slope% row vector
                    
                    om(d3 > SIG*d0 | f3 > f0 + x3.*d0*RHO) = false;             % update outer mask
                    if ~any(om) || M == 0  % are we done extrapolating?
                        break; end
                    
                    x1(om) = x2(om); f1(om) = f2(om); d1(om) = d2(om); % move point 2 to point 1    % row vectors
                    x2(om) = x3(om); f2(om) = f3(om); d2(om) = d3(om); % move point 3 to point 2    % row vectors
                    A(om) = 6*(f1(om)-f2(om))+3*(d2(om)+d1(om)).*(x2(om)-x1(om));% make cubic extrap.% row vectors
                    B(om) = 3*(f2(om)-f1(om))-(2*d1(om)+d2(om)).*(x2(om)-x1(om));
                    x3(om) = x1(om)-d1(om).*(x2(om)-x1(om)).^2./(B(om)+sqrt(B(om).*B(om)-A(om).*d1(om).*(x2(om)-x1(om)))); % num. error possible, ok!
                    
                    % implement fixups: num prob | wrong sign | new point beyond extrapolation limit?
                    issue_mask = ~isreal(x3) | isnan(x3) | isinf(x3) | x3 < 0 | x3 > x2*EXT;
                    x3(om & issue_mask) = x2(om & issue_mask)*EXT;
                    % new point too close to previous point?
                    issue_mask = ~issue_mask & (x3 < x2+INT*(x2-x1));
                    x3(om & issue_mask) = x2(om & issue_mask) + INT*(x2(om & issue_mask) - x1(om & issue_mask));
                end                                                       % end extrapolation

                while 1                    
                    om = keep_going & (abs(d3) > -SIG*d0 | f3 > f0+x3.*d0*RHO);    % update outer mask
                    if ~(any(om) && M > 0)
                        break; end
                    tog = d3 > 0 | f3 > f0+x3.*d0*RHO;           % subinterval toggle
                    x4(om&tog) = x3(om&tog); 
                    f4(om&tog) = f3(om&tog); 
                    d4(om&tog) = d3(om&tog);
                    x2(om&~tog) = x3(om&~tog); 
                    f2(om&~tog) = f3(om&~tog); 
                    d2(om&~tog) = d3(om&~tog);                    
                    tog=f4>f0;                                   % interpolation toggle
                    x3(om&tog) = x2(om&tog) - (0.5*d2(om&tog).*(x4(om&tog)-x2(om&tog)).^2) ./ (f4(om&tog)-f2(om&tog)-d2(om&tog).*(x4(om&tog)-x2(om&tog)));  % quadratic interpolation                
                    A(om&~tog) = 6*(f2(om&~tog)-f4(om&~tog))./(x4(om&~tog)-x2(om&~tog)) + 3*(d4(om&~tog)+d2(om&~tog));                    % cubic interpolation
                    B(om&~tog) = 3*(f4(om&~tog)-f2(om&~tog)) - (2*d2(om&~tog)+d4(om&~tog)).*(x4(om&~tog)-x2(om&~tog));
                    x3(om&~tog) = x2(om&~tog) + (sqrt(B(om&~tog).*B(om&~tog) - A(om&~tog).*d2(om&~tog).*(x4(om&~tog)-x2(om&~tog)).^2)-B(om&~tog))./A(om&~tog);  % num. error possible, ok!
                    numprob = isnan(x3) | isinf(x3);             % if we had a numerical problem then bisect
                    x3(om&numprob) = (x2(om&numprob)+x4(om&numprob))/2;

                    x3(om) = max(min(x3(om), x4(om)-INT*(x4(om)-x2(om))),x2(om)+INT*(x4(om)-x2(om)));  % don't accept too close
                    [f3(om) df3(:,om)] = feval(f, X(:,om) + bsxfun(@times,x3(om),s(:,om)), Wz(:,om), Wu(:,om), ccp(:,om), rho); % masked array op
                    mask = om & f3<F0;                                      % keep best values
                    X0(:,mask) = X(:,mask) + bsxfun(@times,x3(mask),s(:,mask)); 
                    F0(mask) = f3(mask); 
                    dF0(:,mask) = df3(:,mask);
                    M = M - 1; i = i + (len<0);                             % count epochs?!
                    d3(om) = sum(df3(:,om).*s(:,om));                                    % new slope% row vector
                end
                
                % --- we are here ---
                
                okay = abs(d3) < -SIG*d0 & f3 < f0+x3.*d0*RHO;

                upd = okay & keep_going;
                X(:,upd) = X(:,upd) + bsxfun(@times,x3(upd),s(:,upd));
                f0(upd) = f3(upd);
                s(:,upd) = bsxfun(@times,sum(df3(:,upd).*df3(:,upd)) - sum(df0(:,upd).*df3(:,upd)) ./ sum(df0(:,upd).*df0(:,upd)),s(:,upd)) - df3(:,upd);
                df0(:,upd) = df3(:,upd);                                % swap derivatives
                d3(upd) = d0(upd);
                d0(upd) = sum(df0(:,upd).*s(:,upd));
                mask = upd & d0 > 0;                         % otherwise use steepest direction
                s(:,mask) = -df0(:,mask);
                d0(mask) = -sum(s(:,mask).*s(:,mask));                
                x3(upd) = x3(upd) .* min(RATIO, d3(upd)./(d0(upd)-realmin));  % slope ratio but max RATIO                
                ls_failed(upd) = false;                              % this line search did not fail
                
                upd = ~okay & keep_going;
                X(:,upd) = X0(:,upd);
                f0(upd) = F0(upd);
                df0(:,upd) = dF0(:,upd);
                
                keep_going(upd) = keep_going(upd) & ~ls_failed(upd); % line search failed twice in a row
                                                                     % so we give up on these guys
                upd = upd & keep_going;
                if i <= abs(len)
                    s(:,upd) = -df0(:,upd);
                    d0(:,upd) = -sum(s(:,upd).*s(:,upd));                             % try steepest
                    x3(upd) = 1./(1-d0(upd));
                    ls_failed(upd) = true;
                end
            end
        end
        
        function layout = dialog_layout_defaults(self)
            layout = {'SignalProcessing.Resampling.SamplingRate','SignalProcessing.SpectralSelection.FrequencySpecification', ...
                'SignalProcessing.EpochExtraction','','Lambdas','WaveDuration','Stepping','ADMMIterations','CgIterations',...
                'LanczosIterations','NumPatterns','SparsityTraceTradeoff'};
        end
        
    end
end

% (turn off a few editor warnings because some actual implementations are missing in this file)
%#ok<*INUSD,*STOUT,*MANU>


