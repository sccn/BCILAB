% Creates an object from the class eegEpoch
% The class eegEpoch can handle epoched data either in channel or component space.
% The class defines methods for ERP computation, single trial wavelet denoising based on Quian and Garcia (2002), 
% and time-frequency analysis. See examples below.
%
% Author: Alejandro Ojeda, Swartz Center for Computational Neuroscience, UCSD, 28-Aug-2012
% email:  alejandro@sccn.ucsd.edu
% 
% Reference:  Quian Quiroga, R., and Garcia, H., Single-trial event-related potentials with wavelet denoising. 
%               Clinical Neurophysiology 114 (2003) 376???390
%
% Before running the examples below, download example_erp.mat from http://code.google.com/p/mobilab/downloads/
% Then copy the code in a new .m file and press Ctrl+R to uncomment 
% %% Loading the data
% load example_erp.mat;
% % data: event locked time series data (EEG channel data or independent componet), dimensions: number of time points X number of trials
% % time: time vector in seconds, dimensions: number of time points X 1
% % channelLabel: label of the channel/or component, example: 'Pz' for channel Pz, 'IC1' for independent component #1
% % condition: name of the condition, example: 'Target', 'Non target', and so on
% % preStimulusMaxLatency: [time(i) time(j)] interval considered as pre-stimulus period, dimension: two elements vector
% 
% %% Constructing the object
% epochDataObj = eegEpoch('data',data,'time',time,'channelLabel',label,'condition',condition,'preStimulusMaxLatency',preStimulusMaxLatency);
%         
% %% Plotting stacked trials and ERP
% plot(epochDataObj);
%         
% %% Time/Frequency Analysis
% % Log-spaced frequencies are recommended because that way you are doing narrow-band decomposition at low frequencies and broad-band at high frequencies.
% % See the help of scal2frq on the Wavelet Toolbox
% fmin = 2;            % 1Hz 
% fmax = 40;           % 40 Hz
% numFreq = 64         % 64 frequencies
% wname = 'cmor1-1.5'; % complex morlet wavelet, see the help of 'cwt' on the Wavelet Toolbox
% plotFlag = true;     % plot the results
% [coefficients,ersp,itc,frequency,time] = waveletTimeFrequencyAnalysis(epochDataObj,wname,fmin,fmax,numFreq,plotFlag);


classdef eegEpoch < epochObject
    properties
        preStimulusMaxLatency
        erp = [];
    end
    methods
        function obj = eegEpoch(varargin)
            obj@epochObject(varargin{:});
            
            if length(varargin) < 7
                [~,loc] = min(abs(obj.timeStamp));
                obj.preStimulusMaxLatency = [1 loc];
            else
                obj.preStimulusMaxLatency = varargin{7};
            end
        end
        %%
        function erp = get.erp(obj)
            if isempty(obj.erp), 
                disp('Computing wavelet-denoised ERP.');
                obj.computeERP;
            end
            erp = obj.erp;
        end
        %%
        function hFigure = plot(obj,channel)
            dim = obj.mmfObj.Format{1,2};
            if length(dim) == 3, Nch = dim(3);else Nch = 1;end
            if nargin < 2, channel = 1:Nch;end
            Nch = length(channel);
            if Nch > 1
                for it=1:Nch, hFigure = plot(obj,channel(it));end
                return
            end
            [~,loc] = min(abs(obj.timeStamp));
            hFigure = figure('Color',[0.93 0.96 1]);
            subplot(211);imagesc(obj.timeStamp,1:dim(2),obj.data(:,obj.sorting,channel)');
            title(['Trials ' obj.label{channel} '  Condition: ' obj.condition]);
            hold on;plot([1 1]*obj.timeStamp(loc),get(gca,'YLim'),'k-.','LineWidth',2);
            subplot(212);plot(obj.timeStamp,obj.erp(:,channel));
            set(gca,'Xlim',obj.timeStamp([1 end]))
            title(['ERP ' obj.label{channel} '  Condition: ' obj.condition]);
            hold on;plot([1 1]*obj.timeStamp(loc),get(gca,'YLim'),'k-.','LineWidth',2);
            xlabel('Time (sec)')
        end
        %%
        function rmThis = detectOutliers(obj,threshold,plotFlag)
            if nargin < 2, threshold = 0.99;end
            if nargin < 3, plotFlag = false;end
            dim = size(obj.data);
            X = reshape(permute(obj.data,[1 3 2]),[prod(dim([1 3])),dim(2)]);
            D = pdist(X');
            try
                Y = mdscale(D,3);
            catch
                Y = cmdscale(D,3);
            end
            r = sqrt(sum(Y.^2,2));
            th = raylinv(threshold, raylfit(r));
            rmThis = r>th;
            
            if plotFlag
                figure('Color',[0.93 0.96 1]);hold on;
                scatter(Y(:,1),Y(:,2),'.','linewidth',2);
                scatter(Y(rmThis,1),Y(rmThis,2),'r.','linewidth',2);
                title('MDS trials');grid on;
                if any(rmThis), legend({'normal' 'outliers'});
                else legend({'normal'});
                end
                axis xy
            end
        end
        %%
        function removeOutliers(obj,rmThis)
            persistent flag
            if ~isempty(flag), disp('You have removed outliers already.');return;end
            if nargin < 2, rmThis = detectOuliers(obj,0.95);end
            data = obj.data(:,~rmThis,:);
            fid = fopen(obj.binFile,'w');fwrite(fid,data(:),class(data));fclose(fid);
            obj.mmfObj = memmapfile(obj.binFile,'Format',{class(data) size(data) 'x'},'Writable',true);
            flag = 1;
        end
        %%
        function sortingByTrialSimilarity(obj)
            dim = size(obj.data);
            % dataSVD = svdDenoising4ERP(obj.data(:,:,1),8);
            % X = zscore(dataSVD);
            X = zscore(obj.data);
            if length(obj.label) > 1
                X = permute(X,[1 3 2]);
                X = reshape(X,[dim(1)*dim(3) dim(2)]);                
            end
            D= pdist(X');
            Y = mdscale(D,2);
            r = sqrt(sum(Y.^2,2));
            [~,obj.sorting] = sort(r);
        end
        %%
        function erp = computeERP(obj,alpha)
            if nargin < 2, alpha = 0.95;end
            data = waveletDenoising(obj,alpha);
            dim = size(data);
            erp = zeros(dim(1),size(data,3));
            for it=1:size(data,3)
                erp(:,it) = geometric_median( squeeze(data(:,:,it))')';
            end 
            obj.erp = erp;
        end
        %%
        function [fdata, bootstat] = waveletDenoising(obj,alpha,channel)
            if nargin < 2, alpha = 0.95;end
            dim = obj.mmfObj.Format{1,2};
            if length(dim) == 3, Nch = dim(3);else Nch = 1;end
            if nargin < 3, channel = 1:Nch;end
            
            Nch = length(channel);
            if Nch > 1
                fdata = obj.data;
                for it=1:Nch
                    disp(['Channel ' num2str(it)]);
                    fdata(:,:,channel(it)) = waveletDenoising(obj,alpha,channel(it));
                end
                return
            end
            data = obj.data(:,:,channel);
            mu = mean(data(obj.preStimulusMaxLatency(1):obj.preStimulusMaxLatency(2),:));
            X = bsxfun(@minus,data,mu);
            erp = mean(X,2);  %#ok
            
            % scales = logspace(0.59,2.8,64);
            dt = diff(obj.timeStamp([1 2]));
            dim = size(X);
            %--
            s0  = 4*dt;  ds = 0.25;  NbSc = 32;
            wname = 'morl';
            %scales = freq2scales(1, 1./dt/2, 64, wname, dt);
            SIG = {erp,dt}; %#ok
            SCA = {s0,ds,NbSc,'pow',2};
            WAV = {wname,4};
            %--
            
            %-- computing wavelet coefficients
            cwtStruct = cwtft(SIG,'scales',SCA,'wavelet',WAV);             
 
            %-- computing ITC
            cwtStructTmp = cwtft({X,dt},'scales',SCA,'wavelet',WAV);
            nf = size(cwtStructTmp.cfs,1);
            coef = reshape(cwtStructTmp.cfs,[nf dim]);
            P = coef./abs(coef);
            itc = squeeze(abs(mean(P,3)));
            
            coef = permute(coef,[3 1 2]);
            coef = reshape(coef,[dim(2) nf*dim(1)]);
            dim = [dim(2) nf dim(1)];
            nboots = 5;
            bootstat = bootstrp(nboots,@boots_itc,coef,ones(dim(1),1)*dim);
            th = raylinv(alpha, raylfit(bootstat(:)));
            if th < 1
                Iitc = itc(:) > th;
                I = Iitc;
                [x,~] = ind2sub(size(cwtStruct.cfs),find(I));
                ux = unique(x);
                [~,loc]= min(abs(((1:NbSc)'*ones(1,length(ux)) - ones(NbSc,1)*ux')));
                Nl = length(loc);
                nelem = zeros(Nl,1);
                for it=1:Nl, nelem(it) = sum(x==loc(it));end
                rmThis = nelem < 0.25*median(nelem);
                if sum(rmThis) < Nl, ux(rmThis) = [];end
                ux(ux<4) = [];
                I = ismember(1:NbSc,ux);
                I = I(:)*ones(1,length(cwtStruct.omega));
                I = logical(I(:));
            else
                I = true(numel(cwtStruct.cfs),1);
            end
            %--
            
            fdata = X;
            for it=1:dim(1)
                cwtStruct = cwtft({X(:,it),dt},'scales',SCA,'wavelet',WAV);
                cwtStruct.cfs(~I) = 0;
                fdata(:,it) = icwtft(cwtStruct,'signal',{X(:,it),dt})';
            end
        end
        %%
        function [coefficients,ersp,itc,frequency,time] = waveletTimeFrequencyAnalysis(obj,wname,fmin,fmax,numFreq,plotFlag,numberOfBoundarySamples,multCompCorrectionMethod, varargin)
            T = diff(obj(1).timeStamp([1 2]));
            if nargin < 2, wname = 'cmor1-1.5';end
            if nargin < 3, fmin = 2;end
            if nargin < 4, fmax = 1/T/2;end
            if nargin < 5, numFreq = 64;end
            if nargin < 6, plotFlag = true;end
            if nargin < 7, numberOfBoundarySamples = 0;end
            if nargin < 8, multCompCorrectionMethod = 'none';end
            
            Nsubjects = length(obj);
            if Nsubjects > 1
                [ersp,itc,frequency,time] = subjectLevelWaveletTimeFrequencyAnalysis(obj,wname,fmin,fmax,numFreq,plotFlag, multCompCorrectionMethod, varargin);
                coefficients = [];
                return
            end
                        
            data = obj.mmfObj.Data.x;
            dim = size(obj.data);
            data = reshape(data,[size(data,1) prod(dim(2:end))]);
            scales = freq2scales(fmin, fmax, numFreq, wname, T);
            frequency = scal2frq(scales,wname,T);
            frequency = fliplr(frequency);
            
            if ~numberOfBoundarySamples
                toCut = round(0.05*length(obj.timeStamp));
            else
                toCut = numberOfBoundarySamples;
            end
            time = obj.timeStamp(toCut:end-toCut-1);
            %-- computing wavelet coefficients
            coefficients = zeros([length(scales) dim(1) prod(dim(2:end))]);
            hwait = waitbar(0,'Computing cwt...','Color',[0.93 0.96 1]);
            prodDim = prod(dim(2:end));
            for it=1:prodDim
               coefficients(:,:,it) = cwt(data(:,it),scales,wname);
               waitbar(it/prodDim,hwait);
            end
            close(hwait);
            
            % fliping frequency dimension
            coefficients = permute(coefficients,[2 1 3]);
            coefficients = reshape(coefficients,[dim(1) length(scales) dim(2:end)]);
            coefficients = flipdim(coefficients,2);
                        
            if toCut > obj.preStimulusMaxLatency(1), t1 = toCut; else t1 = obj.preStimulusMaxLatency(1);end
            if length(obj.timeStamp)-toCut <= obj.preStimulusMaxLatency(2)
                obj.preStimulusMaxLatency(2) = length(obj.timeStamp)-toCut-t1;
                t2 = obj.preStimulusMaxLatency(2);
            else
                t2 = length(obj.timeStamp)-toCut;
            end
            
            coefficientsDB = 10*log10(abs(coefficients).^2+eps);
            base = mean(coefficientsDB(t1:t2,:,:,:));
            coefficients   = coefficients(toCut:end-toCut-1,:,:,:);
            coefficientsDB = 10*log10(abs(coefficients).^2+eps);
            
            ersp = bsxfun(@minus,coefficientsDB,(base)+eps);
            % ersp = squeeze(mean(ersp,3));
            ersp = permute(ersp,[3 1 2 4]);
            finalDim = size(ersp);
            ersp = geometric_median(ersp(:,:));
            ersp = reshape(ersp,finalDim(2:end));
                        
            itc = coefficients./abs(coefficients);
            itc = squeeze(abs(mean(itc,3)));
            
            Nv = length(varargin);
            switch multCompCorrectionMethod
                case 'none'
                    % disp('Not significance test was computed.');
                case 'bootstrap'
                    if Nv < 1, nboot = 1000; else nboot = varargin{1};end
                    if Nv < 2, alpha = 0.05; else alpha = varargin{2};end
                    
                    % ersp
                    coefficientsDB = permute(coefficientsDB,[3 setdiff(1:ndims(coefficientsDB),3)]);
                    dim = size(coefficientsDB);
                    coefficientsDB = reshape(coefficientsDB,[dim(1) prod(dim(2:end))]);
                    bootstat = bootstrp(nboot,@boots_ersp,coefficientsDB,ones(dim(1),1)*[t1 t2],ones(dim(1),1)*dim);
                    bootstat = reshape(bootstat,[nboot dim(2:end)]);
                    ersp = reshape(ersp,[prod(dim(2:3)) length(obj.label)]);
                    I1 = false(prod(dim(2:3)),length(obj.label));
                    I2 = false(prod(dim(2:3)),length(obj.label));
                    for it=1:length(obj.label)
                        tmp = bootstat(:,:,:,it);
                        tmp = reshape(tmp,[nboot prod(dim(2:3))]);
                        maxmin = prctile(tmp,100*[alpha 1-alpha],2);
                        % th = [min(th(:,1)) max(th(:,2))];
                        th(1) = prctile(maxmin(:,1),100*alpha);
                        th(2) = prctile(maxmin(:,2),100*(1-alpha));
                        I = ersp(:,it) > th(1) & ersp(:,it) < th(2);
                        I1(:,it) = I;
                        ersp(I,it) = 0;
                    end
                    ersp = reshape(ersp,[dim(2:3) length(obj.label)]);
                    
                    % itc
                    coefficientsTmp = permute(coefficients,[3 setdiff(1:ndims(coefficients),3)]);
                    coefficientsTmp = reshape(coefficientsTmp,[dim(1) prod(dim(2:end))]);
                    bootstat = bootstrp(nboot,@boots_itc,coefficientsTmp,ones(dim(1),1)*dim);
                    bootstat = reshape(bootstat,[nboot dim(2:end)]);
                    itc = reshape(itc,[prod(dim(2:3)) length(obj.label)]);
                    
                    for it=1:length(obj.label)
                        % th = raylinv((1-alpha), raylfit(itc(:,it)));
                        tmp = bootstat(:,:,:,it);
                        tmp = reshape(tmp,[nboot prod(dim(2:3))]);
                        th = prctile(tmp,100*(1-alpha),2);
                        th = prctile(th,100*(1-alpha));
                        I = itc(:,it) < th;
                        I2(:,it) = I;
                        itc(I,it) = 0;
                    end
                    itc = reshape(itc,[dim(2:3) length(obj.label)]);
                    
                otherwise
                    error('Unknown method. Stick to bootstrap by now.');
            end
            
            if plotFlag
                G = fspecial('gaussian',[4 4],2);
                ersp_s = ersp;
                itc_s = itc;
                for it=1:length(obj.label)
                    ersp_s(:,:,it) = imfilter(ersp_s(:,:,it),G,'same');
                    itc_s(:,:,it)  = imfilter(itc_s(:,:,it), G,'same');
                    %ersp(I1) = 0;
                    %itc(I2) = 0;
                    
                    eegEpoch.imageLogData(time,frequency,ersp_s(:,:,it));
                    title(['ERSP (dB) ' obj.label{it} '  Condition: ' obj.condition]);
                    
                    strTitle = ['ITC ' obj.label{it} '  Condition: ' obj.condition];
                    eegEpoch.imageLogData(time,frequency,itc_s(:,:,it),strTitle);
                end
            end
            
        end
    end
    methods(Hidden)
        function [t_ersp,t_itc,frequency,time] = subjectLevelWaveletTimeFrequencyAnalysis(obj,wname,fmin,fmax,numFreq,plotFlag, multCompCorrectionMethod, varargin)
            Nsubjects = length(obj);
            if Nsubjects < 2, error('You must input an array of eegEpoch objects, each element in the array containing single subject data.');end
           
            T = diff(obj(1).timeStamp([1 2]));
            if nargin < 2, wname = 'cmor1-1.5';end
            if nargin < 3, fmin = 2;end
            if nargin < 4, fmax = 1/T/2;end
            if nargin < 5, numFreq = 64;end
            if nargin < 6, plotFlag = true;end
            if nargin < 7, multCompCorrectionMethod = 'none';end
            
            [~,ersp,itc,frequency,time] = waveletTimeFrequencyAnalysis(obj(1),wname,fmin,fmax,numFreq,false);
            ersp = repmat(ersp,[1 1 Nsubjects]);
            itc = repmat(itc,[1 1 Nsubjects]);
            for it=2:Nsubjects
                [~,ersp(:,:,it),itc(:,:,it)] = waveletTimeFrequencyAnalysis(obj(it),wname,fmin,fmax,numFreq,false);
                if ~mod(it,10), fprintf(' %i%',round(100*it/Nsubjects));end
            end
            fprintf('\n');
            ersp = permute(ersp,[3 1 2]);
            itc = permute(itc,[3 1 2]);
            t_ersp = tStudent2Dmap(ersp);
            t_itc = tStudent2Dmap(itc);
            Nv = length(varargin);
            
            switch multCompCorrectionMethod
                case 'bootstrap'
                    if Nv < 1, nboot = 1000; else nboot = varargin{1};end
                    if Nv < 2, alpha = 0.05; else alpha = varargin{2};end
                    if Nv < 3, tail = 'both';else tail = varargin{3}; end
                                        
                    % ersp
                    bootstat{1} = bootstrp(nboot,@tStudent2Dmap,ersp);
                    bootstat{2} = bootstrp(nboot,@tStudent2Dmap,itc);
                    
                    switch tail
                        case 'both'
                            th = prctile(bootstat{1},100*[alpha 1-alpha],2);
                            th = [min(th(:,1)) max(th(:,2))];
                            I = t_ersp > th(1) & t_ersp < th(2);
                            t_ersp(I) = 0; 
                                           
                            th = prctile(bootstat{2},100*[alpha 1-alpha],2);
                            th = [min(th(:,1)) max(th(:,2))];
                            I = t_itc > th(1) & t_itc < th(2);
                            t_itc(I) = 0; 
                            pval = 100*[alpha 1-alpha];
                            
                        case 'right'
                            th = prctile(bootstat{1},100*alpha,2);
                            th = max(th);
                            I = t_ersp < th;
                            t_ersp(I) = 0; 
                                           
                            th = prctile(bootstat{2},100*alpha,2);
                            th = max(th);
                            I = t_itc < th;
                            t_itc(I) = 0; 
                            pval = 100*alpha;
                            
                        case 'left'
                            th = prctile(bootstat{1},100*(1-alpha),2);
                            th = min(th);
                            I = t_ersp > th;
                            t_ersp(I) = 0; 
                                           
                            th = prctile(bootstat{2},100*(1-alpha),2);
                            th = min(th);
                            I = t_itc > th;
                            t_itc(I) = 0; 
                            pval = 100*alpha;
                        otherwise
                            error('Wrong tail, select from: ''both'', or ''right''.');
                    end
                otherwise
                    error('Unknown method. Stick to bootstrap by now.');
            end
            if plotFlag
                strTitle = ['T-ERSP, pval = [ ' num2str(pval)  '],  ' obj(1).label{1} '  Condition: ' obj(1).condition];
                imageLogData(time,frequency,t_ersp,strTitle);
                
                strTitle = ['T-ITC, pval = [ ' num2str(pval)  '],  ' obj(1).label{1} '  Condition: ' obj(1).condition];
                imageLogData(time,frequency,t_itc,strTitle);
            end
        end
    end
    methods(Static)
        function imageLogData(time,frequency,data,strTitle)
            if nargin < 4, strTitle = '';end
            figure('Color',[0.93 0.96 1]);
            imagesc(time,log10(frequency),data');
            hAxes = gca;
            tick = get(hAxes,'Ytick');
            fval = 10.^tick;
            Nf = length(tick);
            yLabel = cell(Nf,1);
            fval(fval >= 10) = round(fval(fval >= 10));
            for it=1:Nf, yLabel{it} = num2str(fval(it),3);end
            mx = max(data(:));
            if min(data(:)) < 0,
                mn = -mx;
            else
                mn = min(data(:));
            end
            set(hAxes,'YDir','normal','Ytick',tick,'YTickLabel',yLabel,'CLim',[mn mx]);
            [~,loc] = min(abs(time));
            hold(hAxes,'on');plot([1 1]*time(loc),get(hAxes,'YLim'),'k-.','LineWidth',2);
            xlabel('Time (sec)');
            ylabel('Frequency (Hz)');
            title(strTitle)
            colorbar;
        end
    end
end

%-
function ersp = boots_ersp(coefficientsDB,preStimulusLatency,dim)
coefficientsDB = reshape(coefficientsDB,dim(1,:));
base = mean(coefficientsDB(:,preStimulusLatency(1,1):preStimulusLatency(1,2),:,:),2);
ersp = bsxfun(@minus,coefficientsDB,base+eps);
ersp = squeeze(mean(ersp));
end

%--
function itc = boots_itc(coefficients,dim)
coefficients = reshape(coefficients,dim(1,:));
itc = coefficients./abs(coefficients);
itc = squeeze(abs(mean(itc)));
end

%--
function [geometricMedian,convergenceHistory,weights]= geometric_median(x, varargin)
% geometricMedian = geometric_median(x, {key, value pairs})
% Input
%
%   x   is an N x M matrix, representing N observations of a M-dimensional matrix.
%
% Key, value pairs
%
%   initialGuess        is optional an 1 x M matrix, representing the initial guess for the gemetrix median
%
%   tolerance           an scalar value. It is the maximum relative change in geometricMedian vector (size of the change in
%                       the last iteration divided by the size of the geometricMedian vector) that makes the
%                       algorithm to continue to the next iteration. If relative change is less than tolerance, it is assumed
%                       that convergence is achieved.
%                       had a relative change more than tolerance then more iterations are performed.
%                       default = 1e-4.
%
% Output
%   geometricMedian     is an 1 x m matrix.
%   convergenceHistory  shows the value of maximum relative chage, which is compared to tolerance in
%                       each iteration.

% use mean as the median as an initial guess if none is provided.

inputOptions = finputcheck(varargin, ...
    {'initialGuess'         'real'  [] mean(x);...
    'tolerance'             'real' [0 1] 1e-4;...
    'maxNumberOfIterations' 'integer' [1 Inf]  1000;...
    });

geometricMedian = inputOptions.initialGuess;
for i=1:inputOptions.maxNumberOfIterations
    lastGeometricMedian = geometricMedian;
    differenceToEstimatedMedian = bsxfun(@minus, x, geometricMedian);
    sizeOfDifferenceToEstimatedMedian = (sum(differenceToEstimatedMedian .^2, 2) .^ 0.5);
    oneOverSizeOfDifferenceToEstimatedMedian = 1 ./ sizeOfDifferenceToEstimatedMedian;
    
    % to prevent nans
    oneOverSizeOfDifferenceToEstimatedMedian(isinf(oneOverSizeOfDifferenceToEstimatedMedian)) = 1e20;    
    geometricMedian = sum(bsxfun(@times, x , oneOverSizeOfDifferenceToEstimatedMedian)) / sum(oneOverSizeOfDifferenceToEstimatedMedian);
    %maxRelativeChange = max(max(abs(lastGeometricMedian - geometricMedian)) ./ abs(geometricMedian));
    maxRelativeChange = (sum((lastGeometricMedian - geometricMedian).^2) / sum(geometricMedian.^2)) .^ 0.5;
    
    if nargout > 1, convergenceHistory(i) = maxRelativeChange;end;
    if (maxRelativeChange < inputOptions.tolerance || isnan(maxRelativeChange)), break;end
end

if nargout > 2
    differenceToEstimatedMedian = bsxfun(@minus, x, geometricMedian);
    sizeOfDifferenceToEstimatedMedian = (sum(differenceToEstimatedMedian .^2, 2) .^ 0.5);
    oneOverSizeOfDifferenceToEstimatedMedian = 1 ./ sizeOfDifferenceToEstimatedMedian;
    
    % to prevent nans
    oneOverSizeOfDifferenceToEstimatedMedian(isinf(oneOverSizeOfDifferenceToEstimatedMedian)) = 1e20;
    weights =  oneOverSizeOfDifferenceToEstimatedMedian / sum(oneOverSizeOfDifferenceToEstimatedMedian);
end
end