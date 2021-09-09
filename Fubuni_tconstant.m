%% Fubuni solution

%% Initialization

clear all; close all; clc;

% Excitation frequency and amplitude
fa  = 1.0e5;    % [Hz]
dpa = 1.0e6;    % [Pa]

% Number of Fubini harmonics in the outer loop
NF = 50;

% Number Bessel summands
NB = 50;

% Speed of sound in given medium
c0 = 1500;   % [m/s]

% Number of spacesteps per simulated time and number of simulated periods
Nspacesteps = 10000;
Nperiods    = 5;

% Time t and distance x for given Nperiods 
t = 2.000e-05;
x = Nperiods*c0/fa;

% Shock formation distance xbar and sigma
xbar = 1.5347;
sigma = linspace(0, x/xbar, Nspacesteps);

% Duration of one excitation period
Ta = 1.0/fa;

% Vectors of dimensionless time and pressure
fatau  = [];
pbydpa = [];


%% Calculation

% Loop over sigma
for i = 1:length(sigma)
    
    % calculating dimensionless time
    tau = t - (sigma(i)*xbar/c0);
    tauAll(i) = tau;  % temporary
    omega = 2*pi*fa;
    
    % Loop over the harmonics of the Fubini solution
    % sumFT is the sum over the Fubini terms (outer loop)
    sumFT = 0;
    for n = 1:NF 
        coeff = 2/(n*sigma(i));
        invcoeff = 1/coeff;
        
        % Loop over the Bessel summands
        % sumBT is the sum of the Bessel summands
        sumBT = 0;
        for m = 0:NB
            % Computing the factorials of the Bessel function
            facm = factorial(m);
            facmplusn = factorial(m+n);
            
            % Computing the Bessel summand and the sum of Bessel summands
            BT = (invcoeff^(2*m+n)) * ((-1)^m)/(facm*facmplusn);
            sumBT = sumBT + BT;
        end
        
        FT = coeff * sin(n*omega*tau) * sumBT;
        sumFT = sumFT + FT;
    end
    
    % dimensionless pressure
    pbydpa(i) = dpa*sumFT/dpa; 
    disp(i)
end


%% Plotting

plot(sigma, pbydpa, 'LineWidth', 1.5);
xlabel('sigma'); ylabel('pbydpa');
grid on;