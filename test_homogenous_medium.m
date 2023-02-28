%% Defining Homogenous Propagation medium

% % create the computational grid
close all;
clear;
 % create the computational grid
Nx = 260;           % number of grid points in the x (row) direction
Ny = 824;           % number of grid points in the y (column) direction

dx = 3e-05;        % grid point spacing in the x direction [m]
dy = 3e-05;        % grid point spacing in the y direction [m]

kgrid = kWaveGrid(Nx, dx, Ny, dy);

% Defining Medium Properties 
medium.sound_speed = ones(Nx,Ny)*157;
medium.density = ones(Nx,Ny)*1055;     % [kg/m^3]
medium.alpha_coeff = ones(Nx,Ny)*0.2055;  % [dB/(MHz^y cm)]
medium.BonA =ones(Nx,Ny)*6.11;
medium.alpha_power = 1.05;

% Defining Time Array
t_end = 5e-5;       % [s]
kgrid.makeTime(medium.sound_speed,[],t_end);

%% Define transducer and sensor and simulate 

arc_pos = [20, 650];         % [grid points]  30, 200, 30, 600
radius = 150;                % [grid points] %change to 70 
diameter = 101;              % [grid points] %change to 129
focus_pos = [150, 650];   % [grid points]   % used to be 
source.p_mask = makeArc([Nx, Ny], arc_pos, radius, diameter, focus_pos);
imshow(source.p_mask)
%%
% define a time varying sinusoidal source
% avg_speed_of_sound = mean(mean(parameters{1}));
% avg_density = mean(mean(parameters{2}));

%I = 1250;
I = 3000;
%I = 300;
%I = 1500;

% p = sqrt(I*avg_speed_of_sound*avg_density);
%p = 46e3;



source_freq = 2.5e6;       % [Hz]2500000
source_mag = 62e3;           % [Pa] try 30 MPa - cite 

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
sensor.record = {'p_max', 'p_rms', 'I_avg', 'I'};

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
imagesc(sensor_data.p_rms)
% Create ylabel
title('57 kPa Source at (19.5 mm, 7.4 mm)')
ylabel('y-position [mm]');

% Create xlabel
xlabel('x-position [mm]');
set(gca,'FontSize',36) % Creates an axes and sets its FontSize to 18
 set(gca, 'FontName', 'Arial')

% axis(axes1,'ij');
% hold(axes1,'off');
% Set the remaining axes properties

xticks([160 320 480 640 800])
xticklabels({'5','10','15','20','25'})
yticks([1 80 150 220])
yticklabels({'8','6','4','2'})
c = colorbar;

%%
colorbar('Ticks',[10000,20000,30000,40000,50000, 60000],...
         'TickLabels',{'10','20','30','40','50', '60'})




