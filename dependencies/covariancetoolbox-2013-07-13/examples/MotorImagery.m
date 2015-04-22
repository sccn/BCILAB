% Example of motor imagery classification
% File data/MotorImagery.mat contain preprocessed data from the subject A09
% from the BCI competition IV, dataset IIa.
% Applied preprocessing step are : 
% 
%   - frequency filtering between 8 and 30 Hz (5-th order butterworth)
%   - time based epoching between 3.5 to 5.5 s after the cue
%   - sample covariance estimation
%
% Classification is trained on training data (epoch index 1:288) and
% applied on test data (epoch index 289:576). Results are evaluated in
% terms of classification accuracy.

clear all;
load('data/MotorImagery');

% Data formating
COVtest = data.data(:,:,data.idxTest);
trueYtest  = data.labels(data.idxTest);

COVtrain = data.data(:,:,data.idxTraining);
Ytrain  = data.labels(data.idxTraining);

%% MDM classification - Multiclass
metric_mean = {'euclid','logeuclid','riemann','ld'};
metric_dist = {'euclid','logeuclid','riemann','ld','kullback'};
acc = zeros(length(metric_mean),length(metric_dist));

for i=1:length(metric_mean)
    for j=1:length(metric_dist)
        Ytest = mdm(COVtest,COVtrain,Ytrain,metric_mean{i},metric_dist{j});
        acc(i,j) = 100*mean(Ytest==trueYtest);
    end
end

disp('------------------------------------------------------------------');
disp('Accuracy (%) - Rows : distance metric, Colums : mean metric');
disp('------------------------------------------------------------------');
displaytable(acc',metric_mean,10,{'.1f'},metric_dist)
disp('------------------------------------------------------------------');

%% Discriminant geodesic filtering + MDM classification - Multiclass
metric_mean = {'euclid','logeuclid','riemann','ld'};
metric_dist = {'euclid','logeuclid','riemann','ld','kullback'};
acc = zeros(length(metric_mean),length(metric_dist));

for i=1:length(metric_mean)
    for j=1:length(metric_dist)
        Ytest = fgmdm(COVtest,COVtrain,Ytrain,metric_mean{i},metric_dist{j});
        acc(i,j) = 100*mean(Ytest==trueYtest);
    end
end
disp('------------------------------------------------------------------');
disp('Accuracy (%) - Rows : distance metric, Colums : mean metric');
disp('------------------------------------------------------------------');
displaytable(acc',metric_mean,10,{'.1f'},metric_dist)
disp('------------------------------------------------------------------');

%% MDM classification - Binary case
metric_mean = 'riemann';
metric_dist = 'riemann';
acc = diag(nan(4,1));

for i=1:4
    for j=i+1:4
        % Select the trials
        ixtrain = (Ytrain==i)|(Ytrain==j);
        ixtest = (trueYtest==i)|(trueYtest==j);
        % Classification
        Ytest = mdm(COVtest(:,:,ixtest),COVtrain(:,:,ixtrain),Ytrain(ixtrain),metric_mean,metric_dist);
        % Accuracy
        acc(i,j) = 100*mean(Ytest==trueYtest(ixtest));
    end
end

disp('------------------------------------------------------------------');
disp('Accuracy (%) - Rows/Colums : Couple of classes');
disp('------------------------------------------------------------------');
displaytable(acc'+acc,{'Right Hand','Left Hand','Foot','Tongue'},10,{'.1f'},{'Right Hand','Left Hand','Foot','Tongue'})
disp('------------------------------------------------------------------');

%% Discriminant geodesic filtering + MDM Classification - Binary case
metric_mean = 'riemann';
metric_dist = 'riemann';
acc = diag(nan(4,1));

for i=1:4
    for j=i+1:4
        % Select the trials
        ixtrain = (Ytrain==i)|(Ytrain==j);
        ixtest = (trueYtest==i)|(trueYtest==j);
        % Classification
        Ytest = fgmdm(COVtest(:,:,ixtest),COVtrain(:,:,ixtrain),Ytrain(ixtrain),metric_mean,metric_dist);
        % Accuracy
        acc(i,j) = 100*mean(Ytest==trueYtest(ixtest));
    end
end

disp('------------------------------------------------------------------');
disp('Accuracy (%) - Rows/Colums : Couple of classes');
disp('------------------------------------------------------------------');
displaytable(acc'+acc,{'Right Hand','Left Hand','Foot','Tongue'},10,{'.1f'},{'Right Hand','Left Hand','Foot','Tongue'})
disp('------------------------------------------------------------------');

%% Kmeans usupervised Classification - Binary case
metric_mean = 'riemann';
metric_dist = 'riemann';
acc = diag(nan(4,1));

% for each couple of classes
for i=1:4
    for j=i+1:4
        % Select the trials
        ixtrain = (Ytrain==i)|(Ytrain==j);
        ixtest = (trueYtest==i)|(trueYtest==j);
        % Classification
        Ytest = kmeanscov(COVtest(:,:,ixtest),COVtrain(:,:,ixtrain),2,metric_mean,metric_dist);
        % Find the right labels
        Classes = unique(trueYtest(ixtest));
        truelabels = (trueYtest(ixtest) == Classes(1))+1;
        % Accuracy
        acc(i,j) = 100*mean(Ytest==truelabels);
        if acc(i,j)<50
            acc(i,j) = 100-acc(i,j);
        end
    end
end

disp('------------------------------------------------------------------');
disp('Accuracy (%) - Rows/Colums : Couple of classes');
disp('------------------------------------------------------------------');
displaytable(acc'+acc,{'Right Hand','Left Hand','Foot','Tongue'},10,{'.1f'},{'Right Hand','Left Hand','Foot','Tongue'})
disp('------------------------------------------------------------------');

%% Tangent Space LDA Classification - Binary case
% the riemannian metric
metric_mean = 'riemann';
% update tangent space for the test data - necessary if test data corresponds to
% another session. by default 0.
update = 0;
acc = diag(nan(4,1));

for i=1:4
    for j=i+1:4
        % Select the trials
        ixtrain = (Ytrain==i)|(Ytrain==j);
        ixtest = (trueYtest==i)|(trueYtest==j);
        % Classification
        Ytest = tslda(COVtest(:,:,ixtest),COVtrain(:,:,ixtrain),Ytrain(ixtrain),metric_mean,update);
        % Accuracy
        acc(i,j) = 100*mean(Ytest==trueYtest(ixtest));
    end
end

disp('------------------------------------------------------------------');
disp('Accuracy (%) - Rows/Colums : Couple of classes');
disp('------------------------------------------------------------------');
displaytable(acc'+acc,{'Right Hand','Left Hand','Foot','Tongue'},10,{'.1f'},{'Right Hand','Left Hand','Foot','Tongue'})
disp('------------------------------------------------------------------');