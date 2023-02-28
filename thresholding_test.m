% loading dicom image 

clear
close all
clc

% X is an m-by-n-by-3 array representing a single-frame truecolor RGB image
dicom = dicomread('A0020_SAG_SPINE');


% this step is so that we only work with information we need 
metadata = dicominfo("A0020_SAG_SPINE");

%[dicom_cropped,rect] = imcrop(dicom); 

I_cropped = imcrop(dicom,[220.51 252.51 839.98 275.98]);
I_gray = rgb2gray(I_cropped);
%%
I_gray = rgb2gray(I_cropped);

I = I_gray; 
K = imadjust(I,[0.3 0.6],[]);
% se = strel('line',3,0);
% Io = imopen(K,se);
% Ie = imerode(Io,se);
% Iobr = imreconstruct(Ie,Io);
% imshow(Iobr)
%%
imshow(K)
r = drawrectangle;
mask = createMask(r);
spine_mask = activecontour(K,mask,20000);
close all;
imshow(spine_mask)
in_hematoma = double(createMask(drawfreehand));
close all; 
%%
figure;
imshow(spine_mask)
[x,y] = getpts;
%%
blood_lower = 42; 
dura_upper = y(1); 
dura_lower = 49; %y(2);
csf_upper = 34; %y(3);
csf_lower = 66; %y(4);
pia_upper = y(5);
pia_lower = 75; %y(6);
cord_upper = y(7);
cord_lower = 268; %y(8);
pia2_upper = 243; %y(9);
pia2_lower = 77; %y(10);
csf2_upper = 249;


close all;

%%
soundSpeed=zeros(size(I_gray,1),size(I_gray,2));%m/s
density=zeros(size(I_gray,1),size(I_gray,2));%kg/m^3

% make the entire background blood for all spots I missed 
mesh_size = size(I_gray);
density(1:mesh_size(1), 1:mesh_size(2)) = 1055; 
soundSpeed(1:mesh_size(1), 1:mesh_size(2)) = 1575;
assignment_matrix = zeros(size(I_gray,1),size(I_gray,2));%m/s


for i = 1:size(I_gray,1) % row 
    for j = 1:size(I_gray,2) % col 
        if i < blood_lower && spine_mask(i,j) == 0
            assignment_matrix(i,j) = 1;
        end
        % dura
        if i > dura_upper && i < dura_lower && spine_mask(i,j) == 1
            soundSpeed(i,j) = 1550;
            density(i,j) = 1500;
            assignment_matrix(i,j) = 1;
        end 
        % csf 
        if j  < 520 
            if i > csf_upper && i < 56 && spine_mask(i,j) == 0
                if assignment_matrix(i,j) == 0
                    soundSpeed(i,j) = 1504.5;
                    density(i,j) = 1007;
                    assignment_matrix(i,j) = 1;
                end 
            end 
        else
            if i > csf_upper && i < csf_lower && spine_mask(i,j) == 0
                if assignment_matrix(i,j) == 0
                    soundSpeed(i,j) = 1504.5;
                    density(i,j) = 1007;
                    assignment_matrix(i,j) = 1;
                end
            end 
        end 
        % pia 
        if i > pia_upper && i < pia_lower && spine_mask(i,j) == 1
            soundSpeed(i,j) = 1548; % replace with pia value 
            density(i,j) = 1502; % replace with pia value 
            assignment_matrix(i,j) = 1;
        end 
        % cord 
        if i > cord_upper && i < cord_lower && spine_mask(i,j) == 0
            if assignment_matrix(i,j) == 0
                soundSpeed(i,j) = 1542;
                density(i,j) = 1075;
                assignment_matrix(i,j) = 1;
            end 
        end 
         % hematoma 
        if in_hematoma(i,j) == 1
            soundSpeed(i,j) = 1558.5;
            density(i,j) = 1065;
            assignment_matrix(i,j) = 1;
        end 
        %second pia 
        if i > pia2_upper && spine_mask(i,j) == 1 
            soundSpeed(i,j) = 1548; % replace with pia value 
            density(i,j) = 1502; % replace with pia value             
        end 

        % all black values under 251 is csf 
        if j > 400
            if i > 260 && spine_mask(i,j) == 0 
                soundSpeed(i,j) = 1504.5;
                density(i,j) = 1007;
            end 
        else
            if i > csf2_upper && spine_mask(i,j) == 0
                soundSpeed(i,j) = 1504.5;
                density(i,j) = 1007;
            end 
        end 

    end 
end 


impedance=density.*soundSpeed;%kg/sec*m^2 (Rayl)

imagesc(I_gray)
heatmap(soundSpeed,'GridVisible','off')
colormap jet



