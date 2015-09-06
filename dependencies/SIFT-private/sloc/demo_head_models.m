% eeglab
try
clear all
close all
clc
catch ME
    disp(ME.message) 
clear all
close all
end

%%
% eeglabPath = fileparts(which('eeglab'));
mobilabPath = 'C:\DevelCore\mobilab'; %'/Users/timmullen/Documents/WORK/Toolboxes/mobilab'; %[eeglabPath filesep 'plugins' filesep 'mobilab'];
% addpath(genpath(mobilabPath));

templateFile = [mobilabPath filesep 'data' filesep 'headModelColin27_11997.mat'];
standardMontage = [mobilabPath filesep 'data' filesep 'standard_1020.elc'];
template = load([mobilabPath filesep 'data' filesep 'headModelColin27_11997.mat']);
individualMontage = [mobilabPath filesep 'data' filesep 'personA_chanlocs.sfp'];
surfFile = 'demo_Colin27_11997.mat';
surfData = template.surfData;
save(surfFile,'surfData');             


%% display template head model: Colin27 + aal atlas + 10/20 montage
[elec,label] = readMontage(standardMontage);
hmObj = headModel('surfaces',surfFile,'atlas',template.atlas,'fiducials',template.fiducials,'channelSpace',elec,'label',label);
plotHeadModel(hmObj);


%% warping individual channel space to template
individualSurfaces = 'demo_Colin27_11997.mat';

hmObj = headModel(individualMontage);
plotMontage(hmObj);

aff = hmObj.warpChannelSpace2Template(templateFile,individualSurfaces,'affine');
plotHeadModel(hmObj);


%% warping template to individual channel space
individualSurfaces = 'warped_demo_Colin27_11997.mat';  % name of the output file

hmObj = headModel(individualMontage);
plotMontage(hmObj);

hmObj.warpTemplate2channelSpace(templateFile,individualSurfaces);
plotHeadModel(hmObj);


%% solving the forward problem with OpenMEEG
conductivity = [0.33 0.022 0.33]; % brain and scalp = 0.33 S/m, skull = 0.022 S/m; these conductivies were taken from
                                  % Valdes-Hernandez et al., 2009, Oostendrop TF, 2000; Wendel and Malmivuo, 2006      
normal2surface = true;
hmObj.computeLeadFieldBEM(conductivity,normal2surface);


%% solving the forward problem with NFT
% conductivity = [0.33 0.022 0.33]; % brain and scalp = 0.33 S/m, skull = 0.022 S/m; these conductivies were taken from
%                                   % Valdes-Hernandez et al., 2009, Oostendrop TF, 2000; Wendel and Malmivuo, 2006      
% normal2surface = true;
% hmObj.computeLeadFieldBEM_NFT(conductivity,normal2surface);



%% solving the inverse problem with sLORETA
% K: lead field matrix
% L: Laplaciian operator
% rmIndices: indices to be removed (the Thalamus)
% surfData(3).vertices: source space

brainStructsToKeep = {'Occipital_Sup_L','Occipital_Mid_L','Occipital_Inf_L'}; 
brainStructsToRemove = setdiff_bc(hmObj.atlas.label,brainStructsToKeep); %{'Thalamus_L','Thalamus_R'};
[sourceSpace,K,L,rmIndices] = getSourceSpace4PEB(hmObj,brainStructsToRemove);
load(hmObj.surfaces)
n = size(surfData(3).vertices,1);
J = zeros(n,1);
Jtrue = J;
Jest = J;
ind = setdiff_bc(1:n,rmIndices);

% simulating some Gaussian sources
x0 = [-81.6328 17.9887 93.8088];  % occipital l
x1 = [30.5384 -12.5882 117.1195]; % frontal superior r
x2 = [0.8625 -58.2585 78.9345];   % supra marginal r
x3 = [25.9583 48.3573 35.8547];   % insula l

d = sqrt(sum((surfData(3).vertices(ind,:) - ones(length(ind),1)*x0).^2,2));
Jtrue(ind) = normpdf(d,0,10);
Vtrue = K*Jtrue(ind);

nlambdas = 100;
plotGCV  = true;
threshold = [25 85]; % threshold = [1 99]; threshold = [];  
Jt = inverseSolutionLoreta(Vtrue,K,L,nlambdas,plotGCV,threshold);  % Jt contains the solution for the source potentials within the chosen space, the remaining elements are zero
Jest(ind) = Jt;
Vest = K*Jest(ind);
hmObj.plotOnModel(Jtrue,Vtrue,'True source');
hmObj.plotOnModel(Jest,Vest,'Estimated source');

disp('adding noise');
% adding noise
snr = 5;
vn = std(Vtrue)/snr*randn(length(Vtrue),1);
Vtrue_noise = Vtrue + vn;

Jt = inverseSolutionLoreta(Vtrue_noise,K,L,nlambdas,plotGCV,threshold);
Jest(ind) = Jt;
Vest = K*Jest(ind);
disp('plotting');
hmObj.plotOnModel(Jtrue,Vtrue,'True source (noisy)');
hmObj.plotOnModel(Jest,Vest,'Estimated source (noisy)');

