function predictions = svmperfclassify(x,y,model,options)

clear fun mex_svm_perf_learn;
clear fun mex_svm_perf_classify;

try
if (nargin<4)
    predictions = mex_svm_perf_classify(x,y,model);
else
    predictions = mex_svm_perf_classify(x,y,model,options);
end;
catch
    fprintf(1,'**************************\n');
    lasterror
    fprintf(1,'**************************\n');
    parm_string
    fprintf(1,'**************************\n');
    rethrow(lasterror);
end;

clear fun mex_svm_perf_learn;
clear fun mex_svm_perf_classify;

end
