%% Non linearity in 1D

%% Initialization

close all; clear all; clc; 

% define initial properties
rho0 = 1000;    % initial ambient density [kg/m^3]
c0 = 1500;      % initial speed of sound [m/s]
alpha0 = 0.75;  % absorption coefficient [dB/(MHz^2 cm)]
CFL = 0.15;     % CFL number
ppw = 150;      % points per wavelength

% kgrid
x_size = 3;                % desired domain size [m]
f_max  = 1e+5;             % desired max freq [Hz]
dx = c0/(ppw * f_max);     % grid spacing in x-direction [m]
Nx = round(x_size/dx);     % num of grid points in x-direction (row)
kgrid = kWaveGrid(Nx, dx);

% medium
medium.sound_speed = c0;
medium.density = rho0;
medium.BonA = 19.48;          % non linearity 
medium.alpha_coeff = alpha0;  % absorption coefficient
medium.alpha_power = 1.5;     % absorption power

% time array
t_end = 2e-3;   % [s]
kgrid.t_array = makeTime(kgrid, medium.sound_speed, CFL, t_end);

% source
source.p_mask = zeros(Nx, 1);
source.p_mask(Nx/2,1) = 1;

source_mag  = 5e+6;   % [Pa]
source_freq = 1e+5;   % [Hz]
source.p = source_mag * sin(2*pi*source_freq * kgrid.t_array);

% sensor
sensor.mask = zeros(Nx, 1);
sensor.mask(:) = 1;

%% Simulation

% input_arg = {'PMLInside', false, 'PlotLayout', true};
input_arg = {'PMLInside', false};
sensor_data = kspaceFirstOrder1D(kgrid, medium, source, sensor, input_arg{:});

%% Post Processing

% obtain index for x = [-0.15, 0.15]
x1 = find(abs(kgrid.x_vec - (-0.15)) < 0.0001);
x2 = find(abs(kgrid.x_vec - (0.15)) < 0.0001);
x_index = x1:x2;

% pa :: curve of amplitude decay (t = 8e-4)
x   = kgrid.x_vec;                            % x distance [m]
xsh = 0.1;                                    % shock formation distance [m]
delta_p = source_mag;                         % excitation amplitude [Hz]
pa = (delta_p * pi) ./ (1 + abs(x/xsh));


% plot simulated sensor data :: P(x) for x = [-0.15, 0.15]m, at t = 8e-5;
figure(2)
plot(kgrid.x_vec(x_index), sensor_data(x_index, round((8e-5)/kgrid.dt)))
xlabel('x-position [m]');  ylabel('Signal Amplitude [Pa]');
title('P(x) at  t = 8e-5');

% plot simulated sensor data :: P(x) for x = [-0.15, 0.15]m, at t = 8e-4;
figure(3)
plot(kgrid.x_vec(x_index), sensor_data(x_index, round((8e-4)/kgrid.dt)))
xlabel('x-position [m]'); ylabel('Signal Amplitude [Pa]');
title('P(x) at  t = 8e-4');

% plot simulated sensor data :: P(x) for x = [-1.5, 1.5]m, at t = 8e-5;
figure(4)
plot(kgrid.x_vec, sensor_data(:, round((8e-5)/kgrid.dt)))
xlabel('x-position [m]');  ylabel('Signal Amplitude [Pa]');
title('P(x) at  t = 8e-5');

% plot simulated sensor data :: P(x) for x = [-1.5, 1.5]m, at t = 8e-4;
figure(5)
plot(kgrid.x_vec, sensor_data(:, round((8e-4)/kgrid.dt)))
xlabel('x-position [m]'); ylabel('Signal Amplitude [Pa]');
title('P(x) at  t = 8e-4');
% plot shock amplitude vs distance for x = [-1.5, 1.5]m, at t = 8e-4;
hold on
plot(kgrid.x_vec, pa, 'r')
hold off
ylim([-0.6e+7, 0.6e+7]);

% ####################################################################### %

figure(6)
% plot simulated sensor data :: P(x) for x = [-0.15, 0.15]m, at t = 8e-5;
subplot(2,2,1)
plot(kgrid.x_vec(x_index), sensor_data(x_index, round((8e-5)/kgrid.dt)))
xlabel('x-position [m]');  ylabel('Signal Amplitude [Pa]');
title('P(x) at  t = 8e-5');

% plot simulated sensor data :: P(x) for x = [-0.15, 0.15]m, at t = 8e-4;
subplot(2,2,2)
plot(kgrid.x_vec(x_index), sensor_data(x_index, round((8e-4)/kgrid.dt)))
xlabel('x-position [m]'); ylabel('Signal Amplitude [Pa]');
title('P(x) at  t = 8e-4');

% plot simulated sensor data :: P(x) for x = [-1.5, 1.5]m, at t = 8e-5;
subplot(2,2,3)
plot(kgrid.x_vec, sensor_data(:, round((8e-5)/kgrid.dt)))
xlabel('x-position [m]');  ylabel('Signal Amplitude [Pa]');
title('P(x) at  t = 8e-5');

% plot simulated sensor data :: P(x) for x = [-1.5, 1.5]m, at t = 8e-4;
subplot(2,2,4)
plot(kgrid.x_vec, sensor_data(:, round((8e-4)/kgrid.dt)))
xlabel('x-position [m]'); ylabel('Signal Amplitude [Pa]');
title('P(x) at  t = 8e-4');