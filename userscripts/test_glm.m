traindata = io_loadset('data:/tutorial/imag_movements1/calib/DanielS001R01.dat');
myapproach = {'CSP' 'Prediction',{'MachineLearning',{'Learner',{'glm','ptype','classification','lambdas',2.^(-4:0.5:4)}}}};
[trainloss,lastmodel,laststats] = bci_train('Data',traindata,'Approach',myapproach,'TargetMarkers',{'StimulusCode_2','StimulusCode_3'},'EvaluationMetric','auc')

myapproach = {'CSP' 'Prediction',{'MachineLearning',{'Learner',{'glm','ptype','regression','lambdas',2.^(-4:0.5:4)}}}};
[trainloss,lastmodel,laststats] = bci_train('Data',traindata,'Approach',myapproach,'TargetMarkers',{'StimulusCode_2','StimulusCode_3'},'EvaluationMetric','auc')

%%

X,y,s2,B,pot,tau, opts,G

