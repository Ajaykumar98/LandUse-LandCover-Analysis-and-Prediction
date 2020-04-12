close all;
clc;
clf;
clear all;
%clear;

I1a = imread('l3-nd44o04-102-064-18mar2012-BAND2.tif');
I1b = imread('l3-nd44o04-102-064-18mar2012-BAND3.tif');
I1c = imread('l3-nd44o04-102-064-18mar2012-BAND4.tif');
I1d = imread('l3-nd44o04-102-064-18mar2012-BAND5.tif');


I1a = imresize(I1a,[128 128]);
I1b = imresize(I1b,[128 128]);
I1c = imresize(I1c,[128 128]);
I1d = imresize(I1d,[128 128]);

I1a=double(I1a);
I1b=double(I1b);
I1c=double(I1c);
I1d=double(I1d);

I2a = imread('l3-nd44o04-102-064-13mar2013-band2.tif');
I2b = imread('l3-nd44o04-102-064-13mar2013-band3.tif');
I2c = imread('l3-nd44o04-102-064-13mar2013-band4.tif');
I2d = imread('l3-nd44o04-102-064-13mar2013-band5.tif');


I2a = imresize(I2a,[128 128]);
I2b = imresize(I2b,[128 128]);
I2c = imresize(I2c,[128 128]);
I2d = imresize(I2d,[128 128]);


I2a=double(I2a);
I2b=double(I2b);
I2c=double(I2c);
I2d=double(I2d);


I3a = imread('L3-ND44O04-102-064-01Apr14-BAND2.tif');
I3b = imread('L3-ND44O04-102-064-01Apr14-BAND3.tif');
I3c = imread('L3-ND44O04-102-064-01Apr14-BAND4.tif');
I3d = imread('L3-ND44O04-102-064-01Apr14-BAND5.tif');


I3a = imresize(I3a,[128 128]);
I3b = imresize(I3b,[128 128]);
I3c = imresize(I3c,[128 128]);
I3d = imresize(I3d,[128 128]);


I3a=double(I3a);
I3b=double(I3b);
I3c=double(I3c);
I3d=double(I3d);

I4a = imread('L3-ND44O04-102-064-02Feb16-BAND2.tif');
I4b = imread('L3-ND44O04-102-064-02Feb16-BAND3.tif');
I4c = imread('L3-ND44O04-102-064-02Feb16-BAND4.tif');
I4d = imread('L3-ND44O04-102-064-02Feb16-BAND5.tif');


I4a = imresize(I4a,[128 128]);
I4b = imresize(I4b,[128 128]);
I4c = imresize(I4c,[128 128]);
I4d = imresize(I4d,[128 128]);


I4a=double(I4a);
I4b=double(I4b);
I4c=double(I4c);
I4d=double(I4d);


Ra = zeros(128,128,'uint16');
Rb = zeros(128,128,'uint16');
Rc = zeros(128,128,'uint16');
Rd = zeros(128,128 ,'uint16');


numFeatures = 1;
numResponses = 1;
numHiddenUnits = 200;

layers = [ ...
    sequenceInputLayer(numFeatures)
    lstmLayer(numHiddenUnits)
    fullyConnectedLayer(numResponses)
    regressionLayer];

options = trainingOptions('adam', ...
    'MaxEpochs',100, ...
    'GradientThreshold',1, ...
    'InitialLearnRate',0.05, ...
    'LearnRateSchedule','piecewise', ...
    'LearnRateDropPeriod',125, ...
    'LearnRateDropFactor',0.2, ...
    'Verbose',0, ...
    'Plots','none');

data=[];

for ii = 1:(size(I1a,1))
    for jj = 1:(size(I1a,2))
        data = [I1a(ii,jj),I2a(ii,jj),I3a(ii,jj),I4a(ii,jj)];
        dataTrain = data;
        mu = mean(dataTrain);
        sig = std(dataTrain);
        dataTrainStandardized = (dataTrain - mu) / sig;
        XTrain = dataTrainStandardized(1:end-1);
        YTrain = dataTrainStandardized(2:end);
        net = trainNetwork(XTrain,YTrain,layers,options);
        %net.trainParam.showWindow=0;
        net = predictAndUpdateState(net,XTrain);
        [net,YPred] = predictAndUpdateState(net,YTrain(3));
        YPred = sig*YPred + mu;
        Ra(ii,jj)=YPred;
        net = resetState(net);
  
        
        data = [I1b(ii,jj),I2b(ii,jj),I3b(ii,jj),I4b(ii,jj)];
        dataTrain = data;
        mu = mean(dataTrain);
        sig = std(dataTrain);
        dataTrainStandardized = (dataTrain - mu) / sig;
        XTrain = dataTrainStandardized(1:end-1);
        YTrain = dataTrainStandardized(2:end);
        net = trainNetwork(XTrain,YTrain,layers,options);
        %net.trainParam.showWindow=0;
        net = predictAndUpdateState(net,XTrain);
        [net,YPred] = predictAndUpdateState(net,YTrain(3));
        YPred = sig*YPred + mu;
        Rb(ii,jj)=YPred;
        net = resetState(net);

        
        
        data = [I1c(ii,jj),I2c(ii,jj),I3c(ii,jj),I4c(ii,jj)];
        dataTrain = data;
        mu = mean(dataTrain);
        sig = std(dataTrain);
        dataTrainStandardized = (dataTrain - mu) / sig;
        XTrain = dataTrainStandardized(1:end-1);
        YTrain = dataTrainStandardized(2:end);
        net = trainNetwork(XTrain,YTrain,layers,options);
        %net.trainParam.showWindow=0;
        net = predictAndUpdateState(net,XTrain);
        [net,YPred] = predictAndUpdateState(net,YTrain(3));
        YPred = sig*YPred + mu;
        Rc(ii,jj)=YPred;
        net = resetState(net);

        
        
        
        data = [I1d(ii,jj),I2d(ii,jj),I3d(ii,jj),I4d(ii,jj)];
        dataTrain = data;
        mu = mean(dataTrain);
        sig = std(dataTrain);
        dataTrainStandardized = (dataTrain - mu) / sig;
        XTrain = dataTrainStandardized(1:end-1);
        YTrain = dataTrainStandardized(2:end);
        net = trainNetwork(XTrain,YTrain,layers,options);
        %net.trainParam.showWindow=0;
        net = predictAndUpdateState(net,XTrain);
        [net,YPred] = predictAndUpdateState(net,YTrain(3));
        YPred = sig*YPred + mu;
        Rd(ii,jj)=YPred;
        net = resetState(net);
        fprintf('i and j values are %d and %d \n',ii,jj);

    end
end


imwrite(Ra,'Prediction1.tif');
imwrite(Rb,'Prediction2.tif');
imwrite(Rc,'Prediction3.tif');
imwrite(Rd,'Prediction4.tif');




