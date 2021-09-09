%% Fubuni solution

%% Initialization

% Excitation frequency and amplitude
fa  = 1.0e5;    % [Hz]
dpa = 1.0e6;    % [Pa]

% Number of Fubini harmonics in the outer loop
NF = 20;

% Number Bessel summands
NB = 20;

% Number of time steps per simulated time and number of simulated periods
Ntimesteps = 400;
Nperiods   = 2;

% Position divided by the shock formation distance
% In the current implementation, the Fubini solution is singular at sigma = 0!
sigma = 0.9;    % sigma = x/xsh

% Duration of one excitation period
Ta = 1.0/fa;

% Vectors of dimensionless time and pressure
fatau  = zeros(Ntimesteps, 1);
pbydpa = zeros(Ntimesteps, 1);


%% Calculation

% time loop
for i = 1:Ntimesteps
    % calculating time and dimensionless time
    tau = i*(Ta*Nperiods)/(Ntimesteps-1);
    fatau(i) = fa*tau;
    omegatau = 2*pi*fatau(i);
    
    % Loop over the harmonics of the Fubini solution
    % sumFT is the sum over the Fubini terms (outer loop)
    sumFT = 0;
    for n = 1:NF 
        coeff = 2/(n*sigma);
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
        
        FT = coeff*sin(n*omegatau)*sumBT;
        sumFT = sumFT + FT;
    end
    
    % dimensionless pressure
    pbydpa(i) = dpa*sumFT/dpa; 
end


%% Plotting

plot(fatau, pbydpa, 'LineWidth', 1.5);
xlabel('fatau');  ylabel('pbydpa');
grid on;