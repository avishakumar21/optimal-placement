function parameters  = acoustic_parameters(dicom_image)

%pos = [0 0 840 276];
%mask = rectangle('Position',pos);


imshow(dicom_image);
r = drawrectangle;
mask = createMask(r);
in_hematoma = mask;


% figure;
% imshow(spine_mask)
% [x,y] = getpts;

% blood_lower = 42; 
% dura_upper = y(1); 
% dura_lower = 49; %y(2);
% csf_upper = 34; %y(3);
% csf_lower = 66; %y(4);
% pia_upper = y(5);
% pia_lower = 75; %y(6);
% cord_upper = y(7);
% cord_lower = 268; %y(8);
% pia2_upper = 243; %y(9);
% pia2_lower = 77; %y(10);
% csf2_upper = 249;
%close all;

blood_lower = 42; 
dura_upper = 24.6; 
dura_lower = 49; %y(2);
csf_upper = 34; %y(3);
csf_lower = 66; %y(4);
pia_upper = 49.5;
pia_lower = 75; %y(6);
cord_upper = 58.35;
cord_lower = 268; %y(8);
pia2_upper = 241; %y(9);
pia2_lower = 77; %y(10);
csf2_upper = 249;



% get location values for all the anatomical structures of interest 

I_gray = dicom_image; 
soundSpeed=zeros(size(I_gray,1),size(I_gray,2));%m/s
density=zeros(size(I_gray,1),size(I_gray,2));%kg/m^3
attenCoeff=zeros(size(I_gray,1),size(I_gray,2));%
BonA=zeros(size(I_gray,1),size(I_gray,2));%
atten=zeros(size(I_gray,1),size(I_gray,2));% 
attenpower=1.05;
% make the entire background blood for all spots I missed 
mesh_size = size(I_gray);
density(1:mesh_size(1), 1:mesh_size(2)) = 1055; 
soundSpeed(1:mesh_size(1), 1:mesh_size(2)) = 1575;
attenCoeff(1:mesh_size(1), 1:mesh_size(2)) = 0.2055;%
BonA(1:mesh_size(1), 1:mesh_size(2)) = 6.11;%
atten(1:mesh_size(1), 1:mesh_size(2)) = 6.2;%
assignment_matrix = zeros(size(I_gray,1),size(I_gray,2));%m/s
spine_mask = dicom_image;

for i = 1:size(I_gray,1) % row 
    for j = 1:size(I_gray,2) % col 
        if i < blood_lower && spine_mask(i,j) == 0
            assignment_matrix(i,j) = 1;
        end
        % dura
        if i > dura_upper && i < dura_lower && spine_mask(i,j) == 1
            spine_mask(i,j) = 32
            soundSpeed(i,j) = 1550;
            density(i,j) = 1500;
            attenCoeff(i,j) = 1.1632;
            BonA(i,j) = 6.72;
            atten(i,j) = 33.5;
            assignment_matrix(i,j) = 1;
        end 
        % csf 
        if j  < 520 
            if i > csf_upper && i < 56 && spine_mask(i,j) == 0
                if assignment_matrix(i,j) == 0
                    soundSpeed(i,j) = 1504.5;
                    density(i,j) = 1007;
                    attenCoeff(i,j) = 0.00868;
                    BonA(i,j) = 4.96;
                    atten(i,j) = 0.25;
                    assignment_matrix(i,j) = 1;
                end 
            end 
        else
            if i > csf_upper && i < csf_lower && spine_mask(i,j) == 0
                if assignment_matrix(i,j) == 0
                    soundSpeed(i,j) = 1504.5;
                    density(i,j) = 1007;
                    attenCoeff(i,j) = 0.00868;
                    BonA(i,j) = 4.96;
                    atten(i,j) = 0.25;
                    assignment_matrix(i,j) = 1;
                end
            end 
        end 
        % pia 
        if i > pia_upper && i < pia_lower && spine_mask(i,j) == 1
            soundSpeed(i,j) = 1548; % replace with pia value 
            density(i,j) = 1502; % replace with pia value 
            attenCoeff(i,j) = 1.1632;
            BonA(i,j) = 6.72;
            atten(i,j) = 33.5;
            assignment_matrix(i,j) = 1;
        end 
        % cord 
        if i > cord_upper && i < cord_lower && spine_mask(i,j) == 0
            if assignment_matrix(i,j) == 0
                soundSpeed(i,j) = 1542;
                density(i,j) = 1075;
                attenCoeff(i,j) = 5.98147;
                BonA(i,j) = 6.72;
                atten(i,j) = 172.28;
                assignment_matrix(i,j) = 1;
            end 
        end 
         % hematoma 
        if in_hematoma(i,j) == 1
            soundSpeed(i,j) = 1558.5;
            density(i,j) = 1065;
            attenCoeff(i,j) = 3.09;
            BonA(i,j) = 6.42;
            atten(i,j) = 89;
            assignment_matrix(i,j) = 1;
        end 
        %second pia 
        if i > pia2_upper && spine_mask(i,j) == 1 
            soundSpeed(i,j) = 1548; % replace with pia value 
            density(i,j) = 1502; % replace with pia value   
            attenCoeff(i,j) = 1.1632;
            atten(i,j) = 33.5;
            BonA(i,j) = 6.72;

        end 

        % all black values under 251 is csf 
        if j > 400
            if i > 260 && spine_mask(i,j) == 0 
                soundSpeed(i,j) = 1504.5;
                density(i,j) = 1007;
                attenCoeff(i,j) = 0.00868;
                atten(i,j) = 0.25;
                BonA(i,j) = 4.96;

            end 
        else
            if i > csf2_upper && spine_mask(i,j) == 0
                soundSpeed(i,j) = 1504.5;
                density(i,j) = 1007;
                attenCoeff(i,j) = 0.00868;
                atten(i,j) = 0.25;
                BonA(i,j) = 4.96;

            end 
        end 

    end 
end 


impedance=density.*soundSpeed;%kg/sec*m^2 (Rayl)

imagesc(I_gray)
 
heatmap(soundSpeed,'GridVisible','off')
colormap jet

parameters=cell(1,6);
parameters{1}=soundSpeed;
parameters{2}=density;
parameters{3}=attenCoeff;
parameters{4}=BonA;
parameters{5}=impedance;
parameters{6}=attenpower;
parameters{7}=atten;
