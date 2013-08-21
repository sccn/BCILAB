function model = svmperflearn(x,y,parm_string)

clear fun mex_svm_perf_learn;
clear fun mex_svm_perf_classify;

try
        model = mex_svm_perf_learn(x,y,parm_string);
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