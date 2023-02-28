

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
I_cropped = imcrop(dicom,[220.51 252.51 839.98 275.98]);

I_gray = rgb2gray(I_cropped);

I = I_gray; 
K = imadjust(I,[0.3 0.7],[]);
% se = strel('line',3,0);
% Io = imopen(K,se);
% 
% Ie = imerode(Io,se);
% Iobr = imreconstruct(Ie,Io);
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
% attenCoeff=zeros(size(I_gray,1),size(I_gray,2));
% AonB=zeros(size(I_gray,1),size(I_gray,2)); not sure what this is 
% attenpower=2; % not sure what this is 

% get location values for all the anatomical structures of interest 

soundSpeed=zeros(size(I_gray,1),size(I_gray,2));%m/s
density=zeros(size(I_gray,1),size(I_gray,2));%kg/m^3
% ADD ATTEN COEF AND AONB
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

out=cell(1,3);
out{1}=soundSpeed;
out{2}=density;
% out{3}=attenCoeff;
% out{4}=AonB;
out{3}=impedance;
% out{6}=attenpower;

%% Defining Heterogenous Propagation medium

% make k - grid the size of medium 

[x_len,y_len] = size(I_gray);

% need to make x_len and y_len a multiple of 2 

% create the computational grid
Nx = 276;           % number of grid points in the x (row) direction
Ny = 840;           % number of grid points in the y (column) direction

dx = 0.1e-3;        % grid point spacing in the x direction [m]
dy = 0.1e-3;        % grid point spacing in the y direction [m]

% the size is 0.0414 m (4.14 cm) in the x directions and 0.025 m (2.5 cm) in the y direction 
% NEED TO CHECK THIS IN MICRODICOM 
kgrid = kWaveGrid(Nx, dx, Ny, dy);

%% Defining Medium Properties 

mediumSpace=zeros(Nx,Ny);
medium.sound_speed = ones(Nx,Ny)*1480;
medium.density = ones(Nx,Ny)*1000;     % [kg/m^3]
%medium.alpha_coeff = ones(Nx,Ny)*.002;  % [dB/(MHz^y cm)]
%medium.BonA =ones(Nx,Ny)*5.2;

% what is this step for? 
%medium.sound_speed(:,1)=0;medium.sound_speed(:,213:end)=761;
%medium.density(:,1)=0;medium.density(:,213:end)=1.2;
%medium.alpha_coeff(:,1)=.1;

% define the properties of the propagation medium
medium.sound_speed(:,1:size(I_gray,2)) = out{1};
medium.density(:,1:size(I_gray,2)) = out{2};     % [kg/m^3]
%medium.alpha_coeff(:,1:size(MRI,2)) = out{3};  % [dB/(MHz^y cm)]
%medium.BonA(:,1:size(MRI,2)) = out{4};
%medium.alpha_power = out{6};

%% Defining Time Array

% create time array
t_end = 3e-5;       % [s]
makeTime(kgrid, medium.sound_speed,.05,t_end);


%% Defining Binary Sensory Mask by opposing corners

% create a sensor mask covering the entire computational domain using the
% opposing corners of a rectangle

sensor.mask = [1, 1, Nx, Ny].';

% set the record mode to capture the final wave-field and the statistics at
% each sensor point

sensor.record = {'p','p_final', 'p_max_all', 'p_rms'};

%% Defining Time Varying Pressure Source
%FUS Transducer with a aperature diameter of 35 mm and ROC of 24.5

semicircle_radius = 65; % [grid points]
arc = makeCircle(Nx, Ny, Nx/2, Ny/2, semicircle_radius, pi);

% find total number and indices of the grid points constituting the
% semicircle
arc_indices = find(arc == 1);
Nv = length(arc_indices);


% calculate angles between grid points in the arc and the centre of the
% grid
arc_angles = atan((kgrid.y(arc_indices)) ./ kgrid.x(arc_indices));


% sort the angles into ascending order, and adjust the indices accordingly
[sorted_arc_angles, sorted_index] = sort(arc_angles);
sorted_arc_indices = arc_indices(sorted_index);

% divide the semicircle into Ne separate sensor elements
Ne = 13;
source.p_mask = zeros(Nx, Ny);
for loop = 1:Ne
    
    % get the indices of the grid points belonging to the current element
    % (there is a two grid point gap between the elements)
    voxel_indices = sorted_arc_indices(floor((loop - 1) * Nv / Ne) + ...
        2:floor(loop * Nv / Ne) - 1);
    
    % add the element to the sensor.mask
    sensor.p_mask(voxel_indices) = 1;
    
end
% source.p_mask = zeros(Nx, Ny);
% source.p_mask(1, :) = 1;

%input_args = {'PMLAlpha', [2, 0], 'DisplayMask', display_mask, 'PlotScale', [-0.75, 0.75]};



% radius = round(24.5/(1000*dx));                % [grid points]
% diameter = 2*floor(35/(1000*dx)/2)+1;          % [grid points]
% 
% arclength=2*radius*asin(diameter/(2*radius));
% theta=arclength*180/(pi*radius);
% 
% x2=-radius*sind(theta);
% y2=radius*cosd(theta);
% 
% % Define focal spot               
% focus_pos = [50 82];   % [grid points]
% % define a curved transducer element
% arc_pos = [50 82+radius];        % [grid points]  
% 
% %Determine distance to intended focal point from each element
% DisttoFoc=sqrt((dx*(arc_pos(:,1)-focus_pos(1))).^2+...
%     (dy*(arc_pos(:,2)-focus_pos(2))).^2);
% 
% TimetoFoc=DisttoFoc./medium.sound_speed(1,1);  % [s]
% 
% [~,source.p_mask] = makeMultiArc([Nx, Ny], arc_pos, radius, diameter, focus_pos);

figure
hold on
imagesc(flip(mediumSpace+source.p_mask))
axis tight
title('Simulation Space')
xlabel('mm')
ylabel('mm')

%%
% define a time varying sinusoidal source

intensity=10;%W/cm^2

% source_freq = 515000;       % [Hz]
% source_mag = 100*sqrt(intensity*(medium.sound_speed(1,1)*medium.density(1,1)));           % [Pa]
% source.p = zeros(size(arc_pos,1),length(kgrid.t_array));




% define the initial pressure distribution
source.p0 = zeros(Nx, Ny);
source.p0(39:41, :) = 2;
 
% turn off the PML in the y-direction
input_args = {'PMLAlpha', [2, 0]};


% define a time varying sinusoidal source (instead of an initial pressure)
source_freq = 12e6;     % [Hz]
source_mag = 0.25;      % [Pa]
source.p = source_mag * sin(2 * pi * source_freq * kgrid.t_array);

% define source mask and force to be binary
source.p_mask = source.p0;
source.p_mask(source.p_mask ~= 0) = 1; 

% remove initial pressure field
source = rmfield(source, 'p0');

% run the simulation
sensor_data1 = kspaceFirstOrder2D(kgrid, medium, source, sensor, input_args{:});

% NumPeriods=TimetoFoc/(1/source_freq);
% 
% phase=2*pi*rem(NumPeriods,1);
% 
% source.p(1,:) = source_mag * sin(2 * pi * source_freq * kgrid.t_array - phase(1));

%% Source Filtering

% for i=1:size(arc_pos,1)
% 
% % filter the source to remove high frequencies not supported by the grid
% source.p(i,:) = filterTimeSeries(kgrid, medium, source.p(i,:));
% 
% end



%%
%Running Simulation

sensor.record_start_index=8000;  % what is this 

% create a display mask to display the transducer
display_mask = source.p_mask;

% assign the input options
input_args = {'DisplayMask', display_mask, 'PMLInside', false, 'PlotPML', false};

% run the simulation
sensor_data = kspaceFirstOrder2D(kgrid, medium, source, sensor, input_args{:});

Ifield=.0001*sensor_data.p_max_all.^2./(medium.sound_speed.*medium.density);
figure
pcolor(flip(Ifield))
shading interp
colorbar 
colormap hot


title('Intensity Heatmap (W/cm^2)')
xlabel('mm')
ylabel('mm')

figure
pcolor(flip(sensor_data.p_max_all))
shading interp
colorbar 

title('Pressure Heatmap (Pa)')
xlabel('mm')
ylabel('mm')