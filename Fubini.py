import sys
import matplotlib
import math
import scipy
import numpy as np
import matplotlib.pyplot as plt
from numpy import array
import glob
from collections import namedtuple
import matplotlib.ticker as mticker

# Excitation frequency and amplitude
fa = 1.0e5
dpa = 1.0e6

# Number of Fubini harmonics in the outer loop
NF = 20

# Number Bessel summands
NB = 20

# Number of time steps per simulated time and number of simulated periods
Ntimesteps = 400
Nperiods = 2

# Position divided by the shock formation distance
# In the current implementation, the Fubini solution is singular at sigma = 0!
sigma = 0.9



# Duration of one excitation period
Ta = 1.0/fa

# Vectors of dimensionless time and pressure
fatau = []
pbydpa = []

# Time loop
for i in range(0, Ntimesteps):
  
  # calculating time and dimensionless time
  tau = Ta*float(Nperiods)/float(Ntimesteps-1)*float(i)
  fatau.append(fa*tau)
  omegatau = 2.0*math.pi*fatau[i]
  
  # Loop over the harmonics of the Fubini solution
  # sumFT is the sum over the Fubini terms (outer loop)
  sumFT = 0.0
  for n in range(1, NF+1):
    
    coeff = 2.0/(float(n)*sigma)
    invcoeff = 1.0/coeff
    
    # Loop over the Bessel summands
    # sumBT is the sum of the Bessel summands
    sumBT = 0.0
    for m in range(0, NB+1):
      # Computing the factorials of the Bessel function
      facm = float(math.factorial(m))
      facmplusn = float(math.factorial(m+n))
      
      # Computing the Bessel summand and the sum of Bessel summands
      BT = float((-1)**m)/float(facm*facmplusn)*invcoeff**(2*m+n)
      sumBT = sumBT + BT
    
    FT = coeff*math.sin(float(n)*omegatau)*sumBT
    sumFT = sumFT + FT
      
  pbydpa.append(dpa*sumFT/dpa)

# Plot results
plt.plot(fatau, pbydpa)
plt.show()