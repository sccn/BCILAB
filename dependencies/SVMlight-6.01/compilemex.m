function compilemex( )

try 
    cd src

    fprintf(1,'Compiling mexsvmlearn\n');
    mex -v -O  -DMATLAB_MEX -I. mexsvmlearn.c global.c svm_learn.c svm_common.c svm_hideo.c mexcommon.c mem_clean.c 

    fprintf(1,'Compiling mexsvmclassify\n');
    mex -O  -DMATLAB_MEX -I../kern_search -I. mexsvmclassify.c global.c svm_learn.c svm_common.c svm_hideo.c mexcommon.c mem_clean.c 

    fprintf(1,'Compiling mexsinglekernel\n');
    mex -O  -DMATLAB_MEX -I../kern_search -I. mexsinglekernel.c global.c svm_learn.c svm_common.c svm_hideo.c mexcommon.c mem_clean.c 
 
    fprintf(1,'Compiling mexkernel\n');
    mex -O  -DMATLAB_MEX -I../kern_search -I. mexkernel.c global.c svm_learn.c svm_common.c svm_hideo.c mexcommon.c mem_clean.c 
    
    cd ..
catch
    cd ..
    fprintf(1,'compile failed\n');
end

