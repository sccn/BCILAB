theDir = ['/Users/cdavey/Documents/MyMatlabCode/MFF_Java_Matlab_Release/MatlabExamples'];
%testfilename = 'LLL_02.3_T164_1559.ave';
%testfilename = 'LLL_02.3_T164_1559.seg';
testfilename = 'Long64ChannelWithEvents.mff';

testfilename = [theDir filesep testfilename]


hdr = read_mff_header(testfilename);
 
events = read_mff_event(testfilename);

indType = 'epoch';
startInd = 1;
lastInd = 2;
dataEpoch = read_mff_data(testfilename, indType, startInd, lastInd, []);
figure;
hold on
subplot(2,1,1); 
plot( 1:size(dataEpoch,2), dataEpoch(1,:,1) );
subplot(2,1,2); 
plot( 1:size(dataEpoch,2), dataEpoch(2,:,1) );

indType = 'sample';
startInd = 1;
lastInd = 100;
dataSample = read_mff_data(testfilename, indType, startInd, lastInd, [1 5 10]);
figure;
hold on
subplot(2,1,1); plot( 1:size(dataSample,2), dataSample(1,:) );
subplot(2,1,2); plot( 1:size(dataSample,2), dataSample(2,:) );
