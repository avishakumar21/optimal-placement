% dicom = dicomread('A0020_SAG_SPINE');
% I_cropped = imcrop(dicom,[220.51 252.51 839.98 275.98]);
% 
% I_gray = rgb2gray(I_cropped);
% 
% I = I_gray; 
% K = imadjust(I,[0.3 0.7],[]);
% 
% [rows, columns] = size(K);
% 
% 
% final_im = imageGrid(K);
%createfigure2(imageGrid(I_gray,5,10,15,20,25,2,4,6,8))
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

% thresholding step 
% using values from rect, crop image to only show spine 
I_cropped = imcrop(dicom,[220.51 252.51 839.98 275.98]);

I_gray = rgb2gray(I_cropped);

I = I_gray; 
K = imadjust(I,[0.3 0.7],[]);

imshow(K);
r = drawrectangle;
mask = createMask(r);
spine_mask = activecontour(K,mask,20000);
close all;
spine_mask = double(spine_mask);


%%
imshow(spine_mask)
r = drawrectangle;
mask = createMask(r);
test = activecontour(spine_mask,mask,50, 'edge');
%%
test = double(test);
merge = imfuse(spine_mask, test);
merge = rgb2gray(merge);
imshow(merge)


%%
parameters = acoustic_parameters(spine_mask, test);
[x,y] = size(merge);
sound_speed = ones(x,y);
sound_speed(:,1:y) = parameters{1}(1:x,1:y);

%%

figure;
imagesc(sound_speed);

title('Speed of Sound')
ylabel('x-position [mm]');

% Create xlabel
xlabel('y-position [mm]');
set(gca,'FontSize',36) % Creates an axes and sets its FontSize to 18
% axis(axes1,'ij');
% hold(axes1,'off');
% Set the remaining axes properties

xticks([160 320 480 640 800])
xticklabels({'5','10','15','20','25'})
yticks([1 80 150 220])
yticklabels({'8','6','4','2'})
colorbar;


%%
final_im = imageGrid(merge);




function Igrid = imageGrid(I, varargin)
% imageGrid puts a grid on an image
% Example 1: 
% I = imread('rice.png');
% imageGrid(I);
% Example 2:
% I = imread('cameraman.tif');
% Igrid = imageGrid(I, 'vertLines', 10, 'horzLines', 10, 'lineWidth', 1, 'lineStyle', ':', 'method', 'burn');
% figure
% imshowpair(I, Igrid, 'montage')
% Setting 'method' to 'draw' (default) opens a new figure with the grid
% drawn onto it
% Setting 'method' to 'burn' outputs a new image with the grid burned onto
% it, which requires the computer vision toolbox
% Tested in MATLAB R2018b
% Requires computer vision toolbox
% (c) Clay Swackhamer
% swackhamerclay at gmail.com
% Version 1.2 2018-10-14
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Set defaults
defaultVertLines = 6;
defaultHorzLines = 5;
defaultLineWidth = 1;
defaultLineStyle = '-';
defaultMethod = 'draw';
p = inputParser;
p.CaseSensitive = false;
p.KeepUnmatched = true; %store extra inputs instead of throwing error
addRequired(p, 'I');
addOptional(p, 'vertLines', defaultVertLines);
addOptional(p, 'horzLines', defaultHorzLines);
addOptional(p, 'lineWidth', defaultLineWidth);
addOptional(p, 'lineStyle', defaultLineStyle);
addOptional(p, 'method', defaultMethod);
%Parse inputs
parse(p, I, varargin{:});
vertLines = p.Results.vertLines;
horzLines = p.Results.horzLines;
lineWidth = p.Results.lineWidth;
lineStyle = p.Results.lineStyle;
method = p.Results.method;
drawLogical = strcmpi(method, 'draw'); %True if 'technique' is 'draw'
%Find out where to draw lines
horzTicks = round(linspace(1, size(I,1), horzLines),0); %cols where we will draw vertical lines
vertTicks = round(linspace(1, size(I,2), vertLines),0); %rows where we will draw horizontal lines
maxY = size(I,1); %num pixels in Y
maxX = size(I,2); %num pixels in X
if drawLogical
    %Display the image and draw lines on it
    imshow(I)
    for i = 1:1:vertLines
        line([vertTicks(i), vertTicks(i)], [maxY, 0], 'LineWidth', lineWidth, 'Color', 'g', 'LineStyle', lineStyle); %draw vertical lines
    end
    for i = 1:1:horzLines
        line([0, maxX], [horzTicks(i), horzTicks(i)], 'LineWidth', lineWidth, 'Color', 'g', 'LineStyle', lineStyle); %draw horizontal lines
    end
    Igrid = 'Lines drawn on open figure';
else
    
    %Burn the lines into the image (changes pixels)
    %Put together line coordinates the way that ShapeInserter wants 
    % [x1_line1, y1_line1, x2_line1, y2_line1;...
    %  x1_line2, y1_line2, x2_line2, y2_line2]
    %Vertical lines
    Vlines = zeros(length(vertTicks),4);
    Vlines(:,1) = vertTicks;
    Vlines(:,2) = 0;
    Vlines(:,3) = vertTicks;
    Vlines(:,4) = maxY;
    %Horizontal lines
    Hlines = zeros(length(horzTicks), 4);
    Hlines(:,1) = 0;
    Hlines(:,2) = horzTicks;
    Hlines(:,3) = maxX;
    Hlines(:,4) = horzTicks;
    %Put all lines together
    lines = vertcat(Vlines, Hlines);
    lines = int32(lines);
    %Insert lines and output the modified image
    Igrid = insertShape(I, 'Line', lines, 'LineWidth', lineWidth, 'Color', 'red');
end
ylabel('x-position [mm]');

% Create xlabel
xlabel('y-position [mm]');
set(gca,'FontSize',36) % Creates an axes and sets its FontSize to 18
% axis(axes1,'ij');
% hold(axes1,'off');
% Set the remaining axes properties

xticks([160 320 480 640 800])
xticklabels({'5','10','15','20','25'})
yticks([1 80 150 220])
end
