function sampleNum = mff_nanos2Sample(ns, sampRate)
sampDuration = 1000000000/sampRate;
sampleNum = ns/sampDuration;
sampleNum = fix(sampleNum);

