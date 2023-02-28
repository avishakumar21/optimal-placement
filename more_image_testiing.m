%% loading dicom image 

clear
close all
clc

% X is an m-by-n-by-3 array representing a single-frame truecolor RGB image
dicom = dicomread('A0020_SAG_SPINE');

I_cropped = imcrop(dicom,[224.51 254.51 835.98 301.98]);

I = rgb2gray(I_cropped);

% save I_gray as a jpg 

K = imadjust(I,[0.3 0.7],[]);
se = strel('line',3,0);
Io = imopen(K,se);

Ie = imerode(Io,se);
Iobr = imreconstruct(Ie,Io);
imshow(Iobr)
r = drawrectangle;
mask = createMask(r);
dura = activecontour(Iobr,mask,20000);
%close all
imshow(dura) 
% imshow(Iobr)
% r2 = drawrectangle;
% mask2 = createMask(r2);
% pia = activecontour(Iobr,mask2,3000,'edge');
% close all
% C = imfuse(dura, pia);
% imshow(Iobr)
% r3 = drawrectangle;
% mask3 = createMask(r3);
% hematoma = activecontour(Iobr,mask3,3000,'edge');
% close all
% D = imfuse(C, hematoma);
% imshow(Iobr)
% r4 = drawrectangle;
% mask4 = createMask(r4);
% dura2 = activecontour(Iobr,mask4,3000,'edge');
% close all
% E = imfuse(D, dura2);
% imshow(E)



% imshow(Iobr)
% r5 = drawrectangle;
% mask5 = createMask(r5);
% hematoma = activecontour(Iobr,mask5,300,'edge');
% close all
% F = imfuse(E, hematoma);
% imshow(Iobr)
% r6 = drawrectangle;
% mask6 = createMask(r6);
% dura2 = activecontour(Iobr,mask6,300,'edge');
% close all
% A = imfuse(F, dura2);


