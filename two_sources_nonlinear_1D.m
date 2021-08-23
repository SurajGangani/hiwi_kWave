%% Non linearity in 1D

%% INITIALIZATION

close all; clear all; clc; 

% define initial properties
P_excitation = 30e+6;       % excitation pressure [Pa]
ppw = 200;                  % points per wavelength
CFL = 0.15;                 % CFL number
rho0 = 1000;                % initial ambient(water) density [kg/m^3]
c0 = 1500;                  % initial speed of sound [m/s]
     
% kgrid
x_size = 3;                 % desired domain size [m]
f_max  = 1e+5;              % desired max freq [Hz]
dx = c0/(ppw * f_max);      % grid spacing in x-direction [m]
Nx = round(x_size/dx);      % num of grid points in x-direction (row)
kgrid = kWaveGrid(Nx, dx);

% medium
medium.sound_speed = c0;
medium.density = rho0;
medium.BonA = 5;               % non linearity 
medium.alpha_coeff = 0.75;     % absorption coefficient [dB/(MHz^2 cm)]
medium.alpha_power = 2;        % absorption power

% time array
t_end = 1.000e-3;   % [s]
kgrid.t_array = makeTime(kgrid, medium.sound_speed, CFL, t_end);

% source
source_pos = 0;                  % pos of sources both side of origin [m]
x1_pos = find(abs(kgrid.x_vec - (-source_pos)) < 0.000001);
x2_pos = find(abs(kgrid.x_vec - (+source_pos)) < 0.000001);
source.p_mask = zeros(Nx, 1);
source.p_mask(x1_pos,1) = 1;     % position of 1st source [grid points]
% source.p_mask(x2_pos,1) = 1;   % position of 2nd source [grid points]

source_mag  = P_excitation;   % [Pa]
source_freq = 1e+5;           % [Hz]
source.p(1,:) = source_mag * sin(2*pi*source_freq * kgrid.t_array);   % excitation of 1st source
% source.p(2,:) = source_mag * sin(2*pi*source_freq * kgrid.t_array); % excitation of 2nd source

% sensor
sensor.mask = zeros(Nx, 1);
sensor.mask(:) = 1;
sensor.record = {'p_final'};


%% SIMULATION

input_arg = {'PMLInside', false};
sensor_data = kspaceFirstOrder1D(kgrid, medium, source, sensor, input_arg{:});


%% POST PROCESSING

% obtain index for x = [-0.15, 0.15]
x1 = find(abs(kgrid.x_vec - (-0.15)) < 0.000001);
x2 = find(abs(kgrid.x_vec - (0.15)) < 0.000001);
x_index = x1:x2;

% pa :: curve of amplitude decay (t = 8e-4)
delta_p = source_mag;                        % excitation amplitude [Hz]
% xsh = (0.1*5e+6)/source_mag;               % shock formation distance for Bio Tissue [m]
xsh = (0.3068*5e+6)/source_mag;              % shock formation distance for Water [m]
x   = kgrid.x_vec;                           % x distance [m]
pa  = (delta_p * pi) ./ (1 + abs(x/xsh));


% ==== SPECIAL PLOTS ===================================================================

% 1) plot simulated sensor data :: P(x) for x = [-1.5, 1.5]m
figure(5)
% plot(kgrid.x_vec, sensor_data(:, end))
plot(kgrid.x_vec, sensor_data.p_final)
xlabel('x-position [m]'); ylabel('Signal Amplitude [Pa]');
title( ['P(x)  [' , num2str(ppw) ,' ppw | P = ' , num2str(source_mag/1e6) ,...
    ' MPa | B/A = ' , num2str(medium.BonA) , ' | y = ' , num2str(medium.alpha_power) , ']'] );
% plot shock amplitude vs distance for x = [-1.5, 1.5]m
hold on
plot(x, pa, 'r')  % envelope to wave profile
hold off

% 2) plot normalized simulated sensor data :: P(x) for x = [-1.5, 1.5]m
x_axis = kgrid.x_vec;
% y_axis = sensor_data(:, end)./source_mag; % normalization is applied
y_axis = sensor_data.p_final./source_mag; % normalization is applied
figure(6)
plot(x_axis, y_axis)   
xlabel('x-position [m]'); ylabel('Normalized Signal Amplitude');
title( ['Normalized P(x)  [' , num2str(ppw) ,' ppw | P = ' , num2str(source_mag/1e6) ,...
    ' MPa | B/A = ' , num2str(medium.BonA) , ' | y = ' , num2str(medium.alpha_power) , ']'] );
hold on
plot(x, pa./source_mag, 'r')  % envelope to wave profile
hold off


% === IMPORT VARIABLES INTO CVS ========================================================

% writematrix(x_axis, ['x_axis_' , num2str(source_mag/1e6) , '_MPa.csv']);
% writematrix(y_axis, ['y_axis_' , num2str(source_mag/1e6) , '_MPa.csv']);
