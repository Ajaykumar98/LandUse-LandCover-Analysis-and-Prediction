clc;
clf;
clear;

I1 = imread('Prediction1.tif');
I1s=I1;
I2 =imread('Prediction2.tif');
I3 =imread('Prediction3.tif');
I4 =imread('Prediction4.tif');

F(:,:,1)=I1;
F(:,:,2)=I2;
F(:,:,3)=I3;
F(:,:,4)=I4;
imwrite(F,'finalPrediction.tif');
C=cat(3,I1,I2,I3,I4);
