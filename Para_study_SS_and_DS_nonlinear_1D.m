%% Non linearity in 1D

%% INITIALIZATION

close all; clear all; clc; 

% define initial properties
P_excitation = 5e+6;        % excitation pressure [Pa]
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
medium.BonA = 0;               % non linearity 
medium.alpha_coeff = 0.75;     % absorption coefficient [dB/(MHz^2 cm)]
medium.alpha_power = 2;        % absorption power

% time array
t_end = 1.000e-3;   % [s]
kgrid.t_array = makeTime(kgrid, medium.sound_speed, CFL, t_end);

% shock formation distance
xsh = (0.3068*5e+6)/P_excitation;

% source
source_pos = 0.25*xsh;            % pos of sources both side of origin [m]
x1_pos = find(abs(kgrid.x_vec - (-source_pos)) < 0.00004);
x2_pos = find(abs(kgrid.x_vec - (+source_pos)) < 0.00004);
source.p_mask = zeros(Nx, 1);
% source.p_mask(x1_pos,1) = 1;    % position of 1st source [grid points]
source.p_mask(x2_pos,1) = 1;    % position of 2nd source [grid points]

source_mag  = P_excitation;   % [Pa]
source_freq = 1e+5;           % [Hz]
source.p(1,:) = source_mag * sin(2*pi*source_freq * kgrid.t_array);   % excitation of 1st source
% source.p(2,:) = source_mag * sin(2*pi*source_freq * kgrid.t_array);   % excitation of 2nd source

% sensor
sensor.mask = zeros(Nx, 1);
sensor.mask(round(Nx/2),1) = 1;


%% SIMULATION

input_arg = {'PMLInside', false, 'PlotLayout', true};
sensor_data = kspaceFirstOrder1D(kgrid, medium, source, sensor, input_arg{:});


%% POST PROCESSING

% get maximum pressure
max_pressure = max(sensor_data);
disp(max_pressure)

% ==== SPECIAL PLOTS ===================================================================

[t_sc, scale, prefix] = scaleSI(max(kgrid.t_array(:)));

% 1) plot simulated sensor data :: P(t)
figure(5)
plot(kgrid.t_array*scale, sensor_data)
xlabel(['Time [' prefix 's]']); ylabel('Signal Amplitude [Pa]');
title( ['P(t)  [' , num2str(ppw) ,' ppw | P = ' , num2str(source_mag/1e6) ,...
    ' MPa | B/A = ' , num2str(medium.BonA) , ' | y = ' , num2str(medium.alpha_power) , ']'] );
