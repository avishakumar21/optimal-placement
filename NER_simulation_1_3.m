

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

parameters = acoustic_parameters(spine_mask);


%% Defining Heterogenous Propagation medium

% make k - grid the size of medium 

[x_len,y_len] = size(I_gray);  

% need to make x_len and y_len a multiple of 2 

% x_size = 0.008; %desired domain size in meters 
% c0_min = 1500;
% f_max = 2.5e6;
% points_per_wavelength = 5;
% dx = c0_min/(points_per_wavelength*f_max); % 1.5 e-4 
% Nx_rec = round(x_size/dx);
% 
% y_size = 0.025; %desired domain size in meters 
% dy = c0_min/(points_per_wavelength*f_max); % 1.5 e-4 
% Ny_rec = round(y_size/dy);
% 
% Nx = 54;
% Ny = 162;



% create the computational grid
Nx = 260;           % number of grid points in the x (row) direction
Ny = 824;           % number of grid points in the y (column) direction

% k wave - run on GPU, parallel threads, downsample image 

% Nx = 256; try 248 --> 288 with PML
% Ny = 840; try 824 --> 864 with PML 

dx = 3e-05;        % grid point spacing in the x direction [m]
dy = 3e-05;        % grid point spacing in the y direction [m]

% the size is 8mm in the x direction and 25.13 mm in the y direction
kgrid = kWaveGrid(Nx, dx, Ny, dy);

% Defining Medium Properties 

mediumSpace=zeros(Nx,Ny);
medium.sound_speed = ones(Nx,Ny);
medium.density = ones(Nx,Ny);     % [kg/m^3]
medium.alpha_coeff = ones(Nx,Ny);  % [dB/(MHz^y cm)]
medium.BonA =ones(Nx,Ny);


% define the properties of the propagation medium
medium.sound_speed(:,1:Ny) = parameters{1}(1:Nx,1:Ny);
medium.density(:,1:Ny) = parameters{2}(1:Nx,1:Ny);     % [kg/m^3]
medium.alpha_coeff(:,1:Ny) = parameters{3}(1:Nx,1:Ny);  % [dB/(MHz^y cm)]
medium.BonA(:,1:Ny) = parameters{4}(1:Nx,1:Ny);
medium.alpha_power = parameters{6};


% Defining Time Array

% create time array
t_end = 3e-6;       % [s]
kgrid.makeTime(medium.sound_speed,[],t_end);


%% Define transducer and sensor and simulate 

% define a curved transducer element
arc_pos = [20, 650];         % [grid points]  150 650
radius = 150;                % [grid points] % only change this for focus pos 
% pressure at source mag --> divide A by the surface are 
% check the intensity at focal point 
diameter = 101;              % [grid points] 
focus_pos = [150, 650];      % [grid points]   % this only effects if you have multi element
source.p_mask = makeArc([Nx, Ny], arc_pos, radius, diameter, focus_pos);

cropped_spine = double(spine_mask(1:Nx, 1:Ny));
spine_with_source = imoverlay(cropped_spine,source.p_mask,'white');
%imshow(spine_with_source)
%%
% define a time varying sinusoidal source
avg_speed_of_sound = mean(mean(parameters{1}));
avg_density = mean(mean(parameters{2}));
avg_atten = mean(mean(parameters{7}));

%I = 1250;
I = 3000;
%I = 300;
%I = 1500;

atten_avg = avg_atten; %nepers/meter for blood --> weighted avg  
z = 86*dx;

A = I*exp(-atten_avg*z);


source_p = sqrt(A*avg_speed_of_sound*avg_density);
%source_pressure = source_p/(2 * avg_speed_of_sound * kgrid.dt / dx);
source_pressure = source_p;

%% JOSH RUN THIS SECTION 

source_freq = 2.5e6;       % [Hz]2500000
source_mag = source_pressure;           % [Pa] 
source.p = source_mag * sin(2 * pi * source_freq * kgrid.t_array);

% filter the source to remove any high frequencies not supported by the grid
source.p = filterTimeSeries(kgrid, medium, source.p);
source.p_mode = 'dirichlet';
source.u_mode = 'dirichlet';



%Defining Binary Sensory Mask by opposing corners

% create a sensor mask covering the entire computational domain using the
% opposing corners of a rectangle
sensor.mask = [1, 1, Nx, Ny].';

% set the record mode to capture the final wave-field and the statistics at
% each sensor point
sensor.record = {'p_max', 'p_rms'}; %,'I', 'I_avg'};

% create a display mask to display the transducer
display_mask = source.p_mask;

% assign the input options
input_args = {'DisplayMask', display_mask, 'PlotLayout', true, 'PMLInside', false, 'PlotPML', false};

% run the simulation
sensor_data = kspaceFirstOrder2D(kgrid, medium, source, sensor, input_args{:});

%%

% =========================================================================
% VISUALISATION
% =========================================================================

% plot the simulated sensor data

figure;
normalized = sensor_data.p_max ./ max(sensor_data.p_max(:));
cropped_spine = double(spine_mask(1:Nx, 1:Ny));
spine_with_source = imoverlay(cropped_spine,source.p_mask,'white');
createfigure(spine_with_source,  sensor_data.p_max);

%imagesc(normalized+source.p_mask)



%% proof of concept figure 

%w/cm^2
Ifield = 0.0001.*sensor_data.p_max.^2./(medium.sound_speed.*medium.density);
lower_threshold = 0.03;%from meghana jove submission 
upper_threshold = 30; %100% from High-Intensity Focused Ultrasound Therapy: an Overview for Radiologists
threshold_matrix = zeros(size(sensor_data.p_max));
for i = 1:size(sensor_data.p_max,1) 
    for j = 1:size(sensor_data.p_max,2) 
        if Ifield(i,j) > lower_threshold
            threshold_matrix(i,j) = 5;
        end 
        if Ifield(i,j) > upper_threshold 
            threshold_matrix(i,j) = 10;
        end 

    end 
end 

createfigure(double(spine_mask), threshold_matrix)
title('Intensity (W/cm^2) for 37 kPa Source')
colormap('parula')


