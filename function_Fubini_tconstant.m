%% function of Fubuni solution for constnat time t

%% Initialization

function [sigma, pbydpa] = function_Fubini_tconstant(fa, dpa, c0, t, xbar, x)

    % Number of Fubini harmonics in the outer loop
    NF = 50;

    % Number Bessel summands
    NB = 50;

    % Shock formation distance xbar and sigma
    sigma = x./xbar;

    % Vectors of dimensionless pressure
    % pbydpa = [];
    pbydpa = zeros(1, length(sigma));


    %% Calculation

    % Loop over sigma
    for i = 1:length(sigma)

        % calculating dimensionless time
        tau = t - (sigma(i)*xbar/c0);
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
        if mod(i,100) == 0
            disp( [num2str(i), ' of ', num2str(length(sigma))] )
        end
  
    end


end