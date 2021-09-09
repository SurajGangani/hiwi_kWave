%% Non linearity in 1D

%% INITIALIZATION

close all; clear all; clc; 

% define initial properties
P_excitation = 10e+6;       % excitation pressure [Pa]
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

% shock formation distance (water)
xsh = (0.3068*5e+6)/P_excitation;

% source
source_pos = 0;            % pos of sources both side of origin [m]
x1_pos = find(abs(kgrid.x_vec - (-source_pos)) < 0.000001);
x2_pos = find(abs(kgrid.x_vec - (+source_pos)) < 0.000001);
source.p_mask = zeros(Nx, 1);
% source.p_mask(x1_pos,1) = 1;    % position of 1st source [grid points]
source.p_mask(x2_pos,1) = 1;    % position of 2nd source [grid points]

source_mag  = P_excitation;   % [Pa]
source_freq = 1e+5;           % [Hz]
source.p(1,:) = source_mag * sin(2*pi*source_freq * kgrid.t_array);   % excitation of 1st source
% source.p(2,:) = source_mag * sin(2*pi*source_freq * kgrid.t_array);   % excitation of 2nd source

% get index of grid points at [0.25 0.5 ...]*xsh
num_of_sensor = floor(1.5/xsh);
sensor_pos = zeros(num_of_sensor,1);
j = 1;                               % loop counter
for i = 0:0.25:num_of_sensor
    sensor_pos(j) = find(abs(kgrid.x_vec - (xsh*i)) < 0.00004);
    j = j + 1;
    disp(i)
end

% sensor
% sensor.mask = zeros(Nx, 1);
% sensor.mask(sensor_pos) = 1;
sensor.mask = zeros(Nx, 1);
sensor.mask(:) = 1;
sensor.record = {'p_final', 'p_max'};
sensor.record_start_index = 131000;


%% SIMULATION

input_arg = {'PMLInside', false};
sensor_data = kspaceFirstOrder1D(kgrid, medium, source, sensor, input_arg{:});


%% POST PROCESSING
%% Fay Solution

% pa :: curve of amplitude decay (t = 8e-4)
x   = kgrid.x_vec;                           % x distance [m]
pa  = (source_mag * pi) ./ (1 + abs(x/xsh));

% get maximum pressure
P_max_analytical = pa(sensor_pos);
P_max_simulation = sensor_data.p_max(sensor_pos);

P_ration = P_max_analytical ./ P_max_simulation;


%% Fubini Solution

% Time t to get Fubini Solution
t = 9.8250e-04;

%[sigma,pbydpa] = function_Fubini_tconstant(fa, dpa, c0, t, xbar, x);
[sigma, pbydpa] = function_Fubini_tconstant(source_freq, P_excitation, c0, t, xsh, x);


% Get Fubini solution only between sigma 1 and -1
sigma_fubini_index = find(sigma > -1 & sigma < 1);
sigma_fubini       = sigma(sigma_fubini_index);
pbydpa_fubini      = pbydpa(sigma_fubini_index);

% script to get edge of an wave profile
% Pressure at each peak of Fubini wave profile
pbydpa_fubini_peak = [];
sigma_fubini_peak  = [];
j = 1;
for i = 2:length(pbydpa_fubini)-1
    if pbydpa_fubini(i) > 0
        if (pbydpa_fubini(i-1) < pbydpa_fubini(i)) && (pbydpa_fubini(i) > pbydpa_fubini(i+1))
            pbydpa_fubini_peak(j) = pbydpa_fubini(i);
            sigma_fubini_peak(j)  = sigma_fubini(i);
            j = j + 1;
        end
    end
end


%% Plotting

% 1) plot simulated sensor data :: P(x) for x = [-1.5, 1.5]m
figure(1)
plot(kgrid.x_vec ./ xsh, sensor_data.p_final ./ source_mag)                     % Simulated wave profile
hold on
plot(kgrid.x_vec ./ xsh, sensor_data.p_max ./ source_mag, 'k', 'LineWidth', 1)  % Simulated max pressure
plot(kgrid.x_vec ./ xsh, pa ./ source_mag, 'r', 'LineWidth', 2)                 % Fay solution
% plot(sigma_fubini, pbydpa_fubini, 'LineWidth', 1)                             % Wave profile of Fubini Solution 
plot(sigma_fubini_peak, pbydpa_fubini_peak, 'LineWidth', 2)                     % Max of Fubini Solution 
xline(-3.5, '--k');  xline(-1, '--k');  xline(1, '--k');  xline(3.5, '--k');
hold off
xlabel('sigma = x-position / xsh'); ylabel('Pmax / Pexcitation');
title( ['P   [' , num2str(ppw) ,' ppw | P = ' , num2str(source_mag/1e6) ,...
    ' MPa | B/A = ' , num2str(medium.BonA) , ' | y = ' , num2str(medium.alpha_power) , ']'] );
legend('P final', 'P max', 'Fay Solution', 'Fubini Solution');

