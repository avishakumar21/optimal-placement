%% loading dicom image 

clear
close all
clc

% X is an m-by-n-by-3 array representing a single-frame truecolor RGB image
dicom = dicomread('A0020_SAG_SPINE');

% this step is so that we only work with information we need 
metadata = dicominfo("A0020_SAG_SPINE");

% crop dicom to display only spine anatomy 
% this step for the first time only to get rect vals
%[dicom_cropped,rect] = imcrop(dicom); 

%% thresholding step 
% using values from rect, crop image to only show spine 
I_cropped = imcrop(dicom,[224.51 254.51 835.98 301.98]);

I_gray = rgb2gray(I_cropped);

% save I_gray as a jpg 


% K = medfilt2(I_gray);
% imshow(K)






%numColors = 5;
%L = imsegkmeans(I_gray,numColors);
%B = labeloverlay(I_gray,L);

% [~,threshold] = edge(I_gray,'sobel');
% fudgeFactor = 0.40;
% BWs = edge(I_gray,'sobel',threshold * fudgeFactor);
% imshow(BWs)
% se90 = strel('line',2,90);
% se0 = strel('line',2,0);
% 
% BWsdil = imdilate(I_gray,[se90 se0]);
% BWdfill = imfill(BWsdil,'holes');
%imshow(BWdfill)
% seD = strel('diamond',1);
% BWfinal = imerode(BWdfill,seD);
% BWfinal = imerode(BWfinal,seD);
% %imshow(BWfinal)
% se = strel('disk',3);
% afterOpening = imopen(BWfinal,se);
% %imshow(afterOpening,[]);
% se2 = strel('disk',3);
% closeBW = imclose(afterOpening,se2);
% %figure, imshow(closeBW)

% mask = zeros(size(closeBW));
% mask(5:end-5,5:end-5) = 1;
% bw = activecontour(closeBW,mask,3000);
%imshow(bw)




% imshow(labeloverlay(I_gray,afterOpening))

% se = strel('',5);
% afterOpening = imopen(I_gray,se);
% figure
% imshow(afterOpening,[]);





