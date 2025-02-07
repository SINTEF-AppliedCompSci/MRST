#!/usr/bin/env python3

"""
Functions to generate a solubility table for the H₂-H₂O mixture.

######################################################################
This script employs the RK-EoS (Redlich-Kwong Equation of State) model 
from the paper "Phase behavior and black-oil simulations of hydrogen 
storage in saline aquifers" (Ahmed, E., et al., 2024).
######################################################################
The original EoS model is derived from Raad et al. (2022):
https://doi.org/10.1016/j.ijhydene.2022.09.222
######################################################################
This script extends the thermodynamic models used in the CO₂-water 
mixture (from https://github.com/Simulation-Benchmarks/11thSPE-CSP/tree/main/thermodynamics) 
to the H2-brine system.
"""

from io import StringIO
import argparse
import math
import urllib
import requests  # pylint: disable=import-error
import numpy as np  # pylint: disable=import-error

try:
    import tqdm
    def withProgress(iterable):
        return tqdm.tqdm(iterable)
except ImportError:
    def withProgress(iterable):
        return iterable


class ParameterRange:
    def __init__(self, min: float, max: float, numSamples: int) -> None:
        assert max > min
        assert numSamples >= 2
        self._min = min
        self._max = max
        self._numSamples = numSamples

    def __getitem__(self, i: int) -> float:
        return self._min + i*self.step

    @property
    def min(self) -> float:
        return self._min

    @property
    def max(self) -> float:
        return self._max

    @property
    def numSamples(self) -> float:
        return self._numSamples

    @property
    def step(self) -> float:
        return (self._max - self._min)/float(self._numSamples - 1)


def getH2Densities(temperatures: ParameterRange, pressures: ParameterRange) -> np.ndarray:
    """Return the H2 densities for the given temperature and pressure ranges"""
    print("Retrieving H2 densities for the requested p/T ranges")
    result = np.ndarray(shape=(temperatures.numSamples, pressures.numSamples))
    for i in withProgress(range(temperatures.numSamples)):
        T = temperatures[i]
        query = {
            "Action": "Data",
            "Wide": "on",
            "ID": "C1333740",
            "Type": "IsoTherm",
            "Digits": "12",
            "PLow": f"{pressures.min}",
            "PHigh": f"{pressures.max}",
            "PInc": f"{pressures.step}",
            "T": str(T),
            "RefState": "DEF",
            "TUnit": "C",
            "PUnit": "Pa",
            "DUnit": "kg/m3",
        }
        response = requests.get(
            "https://webbook.nist.gov/cgi/fluid.cgi?" + urllib.parse.urlencode(query)
        )
        response.encoding = "utf-8"
        text = response.text
        values = np.genfromtxt(StringIO(text), delimiter="\t", names=True)
        phaseColIdx = list(values.dtype.names).index("Phase")
        phase = np.genfromtxt(StringIO(text), delimiter="\t", dtype=str, skip_header=1, usecols=[phaseColIdx])
        phaseTransitionForward = np.append((phase[:-1] != phase[1:]), False)
        phaseTransitionBackward = np.insert((phase[1:] != phase[:-1]), 0, False)
        isNotPhaseTransition = np.invert(np.bitwise_or(phaseTransitionForward, phaseTransitionBackward))
        result[i] = values["Density_kgm3"][isNotPhaseTransition]
    return result

def equilibriumConstantH2(temperature):
    """
    Calculate the equilibrium constant (K) of pure hydrogen (H2) at the reference pressure
    based on the temperature, using an empirical correlation from Raad et al. (2022).
    
    Parameters:
    -----------
    temperature : float
        Temperature in Kelvin (K).
        
    Returns:
    --------
    K_H2 : float
        The equilibrium constant for hydrogen (H2) at the given temperature.
    """
    
    # Convert temperature from Kelvin to Celsius
    temperature_C = temperature - 273.15
    
    # Coefficients from Raad et al. (2022) for H2 equilibrium constant
    coefficients_H2 = [2.9947, 4.81373e-3, -5.1773e-5, 1.19052e-7]
    
    # Compute the logarithm of the equilibrium constant for H2 (log K_H2)
    logK_H2 = (coefficients_H2[0] 
               + coefficients_H2[1] * temperature_C 
               + coefficients_H2[2] * temperature_C**2 
               + coefficients_H2[3] * temperature_C**3)
    
    # Return the equilibrium constant (K_H2) by converting from log base 10
    return math.pow(10, logK_H2)

def equilibriumConstantH2O(temperature):
    """
    Calculate the equilibrium constant (K) of pure water (H2O) at the reference pressure
    based on the temperature, using an empirical correlation from Raad et al. (2022).
    
    Parameters:
    -----------
    temperature : float
        Temperature in Kelvin (K).
        
    Returns:
    --------
    K_H2O : float
        The equilibrium constant for water (H2O) at the given temperature.
    """
    
    # Convert temperature from Kelvin to Celsius
    temperature_C = temperature - 273.15
    
    # Coefficients from Raad et al. (2022) for H2O equilibrium constant
    coefficients_H2O = [-2.1817, 2.98e-2, -1.098e-4, 2.048e-7]
    
    # Compute the logarithm of the equilibrium constant for H2O (log K_H2O)
    logK_H2O = (coefficients_H2O[0] 
                + coefficients_H2O[1] * temperature_C 
                + coefficients_H2O[2] * temperature_C**2 
                + coefficients_H2O[3] * temperature_C**3)
    
    # Return the equilibrium constant (K_H2O) by converting from log base 10
    return math.pow(10, logK_H2O)
 

def fugacityCoefficientH2(T, p, rhoH2):
    """
    Calculates the fugacity coefficient (phi) for hydrogen (H2) based on the 
    simplified Redlich-Kwong equation of state (RK-EoS), as described by 
    Raad et al. (2022) and simplified by Ahmed et al. (2024). 

    In this implementation, the water mole fraction in the CO2-rich phase is 
    neglected using simplified mixing rules, which significantly reduces 
    computational complexity while maintaining accuracy in relevant scenarios.

    Parameters:
    -----------
    T : float
        Temperature in Kelvin.
    p : float
        Pressure in Pascals.
    rhoH2 : float
        Density of hydrogen in kg/m^3.

    Returns:
    --------
    phiH2 : float
        Fugacity coefficient for hydrogen (H2).
    """
    # Molar mass of H2 in kg/mol
    molarMassH2 = 2.016e-3

    # Calculate molar volume [cm^3/mol]
    V = 1 / (rhoH2 / molarMassH2) * 1e6  # Convert from m^3/mol to cm^3/mol

    # Convert pressure to bar
    p_bar = p / 1e5  # 1 Pascal = 1e-5 bar

    # Redlich-Kwong equation parameters for hydrogen (H2)
    a_H2 = 1441753.379  # Parameter 'a' for H2 [bar.cm^6/mol^2]
    b_H2 = 18.417       # Parameter 'b' for H2 [cm^3/mol]

    # Universal gas constant [bar.cm^3/mol.K]
    R = 83.1446261815324

    # Calculation of the natural logarithm of the fugacity coefficient ln(Phi_H2)
    lnPhiH2 = (
        math.log(V / (V - b_H2))                  # Logarithmic term for the molar volume
        + b_H2 / (V - b_H2)                       # Term involving molar volume and b_H2
        - 2 * a_H2 / (R * math.pow(T, 1.5) * b_H2) * math.log((V + b_H2) / V)
        + (a_H2 * b_H2) / (R * math.pow(T, 1.5) * b_H2**2) * (
            math.log((V + b_H2) / V) - b_H2 / (V + b_H2)
        )
        - math.log(p_bar * V / (R * T))           # Pressure and volume adjustment term
    )

    # Return the fugacity coefficient by exponentiating ln(Phi_H2)
    return math.exp(lnPhiH2)

def heaviside(x):
    return 0.5 * (np.sign(x) + 1)
    

def activitySalt(saltMolality, temperature, pressure):
    """
    Computes the activity coefficient of salt in a saline solution using 
    an empirical correlation based on temperature, pressure, and salt molality.
    See (Ahmed, E., et al., 2024) for more details

    Parameters:
    -----------
    saltMolality : float
        Molality of salt in the solution [mol/kg].
    temperature : float
        Temperature of the system [K].
    pressure : float
        Pressure of the system [bar].

    Returns:
    --------
    gamma : float
        Activity coefficient of salt in the solution.
    """
    
    # Normalize pressure and temperature to critical values of water.
    normalizedPressure = pressure / 220.064  # Critical pressure of water [bar]
    normalizedTemp = temperature / 647.096   # Critical temperature of water [K]

    # Empirical coefficients
    a1, a2, a3, a4 = 6.5934, 0.0571, 0.3375, 1.2229
    a5, a6, a7, a8 = 0.3448, -2.9756, -0.0725, -4.2276
    a9, a10, a11, a12 = -0.8155, 0.6073, 0.2889, 0.0215

    # Compute natural log of activity coefficient (lnGamma)
    lnGamma = ((a1 + a2 * normalizedPressure / normalizedTemp + a3 / normalizedTemp + a4 * math.log(normalizedTemp)) 
               * math.pow(saltMolality, a5) 
               + (a6 + a7 * normalizedPressure / normalizedTemp + a8 * normalizedTemp + a9 / normalizedTemp) 
               * math.pow(saltMolality, a11) 
               + 2 * a10 * (heaviside(saltMolality) - 0.5) 
               + a12)

    # Return the activity coefficient by exponentiating lnGamma
    return math.exp(lnGamma)


def fugacityCoefficientH2O(T, p, rhoH2):
    """
    Calculates the fugacity coefficient (phi) for water (H2O) based on the 
    simplified Redlich-Kwong equation of state (RK-EoS), as described by 
    Raad et al. (2022) and simplified by Ahmed et al. (2024).

    In this implementation, the water mole fraction in the CO2-rich phase is 
    neglected using simplified mixing rules, allowing for a more efficient 
    calculation while retaining necessary accuracy.

    Parameters:
    -----------
    T : float
        Temperature in Kelvin.
    p : float
        Pressure in Pascals.
    rhoH2 : float
        Density of hydrogen in kg/m^3.

    Returns:
    --------
    phiH2O : float
        Fugacity coefficient for water (H2O).
    """
    # Molar mass of H2 in kg/mol
    molarMassH2 = 2.016e-3

    # Calculate molar volume [cm^3/mol]
    V = 1 / (rhoH2 / molarMassH2) * 1e6  # Convert from m^3/mol to cm^3/mol

    # Convert pressure to bar
    p_bar = p / 1e5  # 1 Pascal = 1e-5 bar

    # Redlich-Kwong equation parameters for H2 and H2O
    a_H2 = 1441753.379    # Parameter 'a' for hydrogen (H2) [bar.cm^6/mol^2]
    a_H2O = 142666655.8   # Parameter 'a' for water (H2O) [bar.cm^6/mol^2]
    a_H2_H2O = math.pow(a_H2 * a_H2O, 0.5)  # Mixed parameter 'a' for H2-H2O system
    b_H2 = 18.417         # Parameter 'b' for hydrogen (H2) [cm^3/mol]
    b_H2O = 21.127        # Parameter 'b' for water (H2O) [cm^3/mol]

    # Universal gas constant [bar.cm^3/mol.K]
    R = 83.1446261815324

    # Calculation of the natural logarithm of the fugacity coefficient ln(Phi_H2O)
    lnPhiH2O = (
        math.log(V / (V - b_H2))                        # Logarithmic term for molar volume
        + b_H2O / (V - b_H2)                           # Term involving molar volume and b_H2
        - 2 * a_H2_H2O / (R * math.pow(T, 1.5) * b_H2) * math.log((V + b_H2) / V)
        + (a_H2 * b_H2O) / (R * math.pow(T, 1.5) * b_H2 * b_H2) * (
            math.log((V + b_H2) / V) - b_H2 / (V + b_H2)
        )
        - math.log(p_bar * V / (R * T))               # Pressure and volume adjustment term
    )

    # Return the fugacity coefficient by exponentiating ln(Phi_H2O)
    return math.exp(lnPhiH2O)


def computeA(T, p, rhoH2):
    """
    Computes the adjusted equilibrium constant for water (H2O) in the presence 
    of hydrogen (H2) based on the given temperature, pressure, and hydrogen density.

    The function applies the following relationships:
    - Adjusts the pressure from a reference of 1 bar to the current pressure.
    - Uses the average partial molar volume of water.
    - Calculates the fugacity coefficient of water in the H2-water system.

    Parameters:
    -----------
    T : float
        Temperature in Kelvin.
    p : float
        Pressure in Pascals.
    rhoH2 : float
        Density of hydrogen (H2) in kg/m^3.

    Returns:
    --------
    A : float
        The adjusted equilibrium constant for the H2O-H2 system.
    """
    # Calculate the pressure difference from the reference pressure (1 bar)
    deltaP = p / 1e5 - 1  # Convert pressure to bar and subtract reference

    # Average partial molar volume of water (H2O) in cm^3/mol
    v_av_H2O = 18.1  
    
    # Universal gas constant [bar.cm^3/(mol.K)]
    R = 83.1446261815324  
    
    # Calculate the equilibrium constant for water (H2O) at the given temperature
    k0_H2O = equilibriumConstantH2O(T)  
    
    # Calculate the fugacity coefficient of water for the H2-water system
    phi_H2O = fugacityCoefficientH2O(T, p, rhoH2)  
    
    # Convert pressure to bar for further calculations
    p_bar = p / 1e5  
    
    # Compute the adjusted equilibrium constant for the H2O-H2 system
    A = k0_H2O / (phi_H2O * p_bar) * math.exp(deltaP * v_av_H2O / (R * T))
    
    return A


def computeB(T, p, rhoH2):
    """
    Computes the coefficient B for the solubility of hydrogen (H2) in water (H2O).

    Parameters:
    -----------
    T : float
        Temperature in Kelvin.
    p : float
        Total pressure in Pascals.
    rhoH2 : float
        Density of hydrogen (H2) in kg/m^3.

    Returns:
    --------
    B : float
        The solubility coefficient B for H2 in H2O.
    """
    # Calculate the deviation in pressure from the reference pressure (1 bar)
    deltaP = p / 1e5 - 1  # Convert pressure to bar and subtract reference
    
    # Average partial molar volume of hydrogen (H2) in cm^3/mol
    v_av_H2 = 26.7  # [cm^3/mol]
    
    # Universal gas constant [bar.cm^3/(mol.K)]
    R = 83.1446261815324  
    
    # Calculate the equilibrium constant for hydrogen (H2) at the given temperature
    k0_H2 = equilibriumConstantH2(T)  
    
    # Calculate the fugacity coefficient of hydrogen for the water-H2 system
    phi_H2 = fugacityCoefficientH2(T, p, rhoH2)  
    
    # Convert total pressure to bar for further calculations
    p_bar = p / 1e5  
    
    # Compute the solubility coefficient B for H2 in H2O
    B = (phi_H2 * p_bar) / (55.508 * k0_H2) * math.exp(-(deltaP * v_av_H2) / (R * T))
    
    return B


def computeBsalt(T, p, rhoH2, ms):
    deltaP = p/1e5 - 1 # pressure range [bar] from p0 = 1 bar to p
    v_av_H2 = 26.7 # average partial molar volume of H2 [cm3/mol]
    R = 83.1446261815324 # universal gas constant [bar.cm3/mol.K]
    k0_H2 = equilibriumConstantH2(T) # equilibrium constant for H2 at 1 bar
    # ms is the salt molality
    gammSalt = activitySalt(ms, T, p/1e5) # salt activity
    phi_H2 = fugacityCoefficientH2(T, p, rhoH2) # fugacity coefficient of H2 for the water-H2 system
    p_bar = p/1e5
    return phi_H2*p_bar/(55.508*k0_H2*gammSalt)*math.exp(-(deltaP*v_av_H2)/(R*T))

parser = argparse.ArgumentParser(
    description="This script generates tables for H2 - H2O solubilities \n"
    "according to Spycher, Pruess, Ennis-King (2003).\n"
)
parser.add_argument(
    "-t1", "--min_temp", required=True, type=float, help="The minimum temperature in degree Celcius."
)
parser.add_argument(
    "-t2", "--max_temp", required=True, type=float, help="The maximum temperature in degree Celcius."
)
parser.add_argument(
    "-nt",
    "--n_temp",
    required=True,
    type=int,
    help="The number of temperature sampling points."
    " min_temp is the first sampling point, max_temp the last.",
)
parser.add_argument(
    "-p1", "--min_press", required=True, type=float, help="The minimum phase pressure in Pascal."
)
parser.add_argument(
    "-p2", "--max_press", required=True, type=float, help="The maximum phase pressure in Pascal."
)
parser.add_argument(
    "-np",
    "--n_press",
    required=True,
    type=int,
    help="The number of pressure sampling points."
    " min_press is the first sampling point, max_press the last.",
)
parser.add_argument(
    "-sm", "--salt_mole", required=True, type=float, help="The salt molality."
)
parser.add_argument(
    "-o", "--output", default="", help="Output filename (default: generated filename)."
)
cmdArgs = vars(parser.parse_args())

minTemp = cmdArgs["min_temp"]
maxTemp = cmdArgs["max_temp"]
minPress = cmdArgs["min_press"]
maxPress = cmdArgs["max_press"]
saltMole = cmdArgs["salt_mole"]
# print pressure values to check if they are correct
print(f"Minimum Pressure (Pa): {minPress}")
print(f"Maximum Pressure (Pa): {maxPress}")
# print temperature values to check if they are correct
print(f"Minimum Temperature (°C): {minTemp}")
print(f"Maximum Temperature (°C): {maxTemp}")
# print salt molality values to check if they are correct
print(f"Salt Molality (mole/Kg): {saltMole}")
pressures = ParameterRange(
    min=cmdArgs["min_press"],
    max=cmdArgs["max_press"],
    numSamples=cmdArgs["n_press"]
)
temperatures = ParameterRange(
    min=cmdArgs["min_temp"],
    max=cmdArgs["max_temp"],
    numSamples=cmdArgs["n_temp"]
)

densities = getH2Densities(temperatures, pressures);

fileName = "solubilities.csv"
# Set output filename, using the default if -o is not set or empty
if cmdArgs["output"]:
    fileName = cmdArgs["output"]
else:
    fileName = f"SolubilityValues_{minPress}_to_{maxPress}_bar_{minTemp:.1f}_to_{maxTemp:.1f}_C_{saltMole:.1f}.csv"
outFile = open(fileName, "w")
outFile.write(f"# This autogenerated file contains solubilities of H2 and H2O in a respective fluid system.\n")
outFile.write("# The values have been calculated by means of (11)-(14) in https://doi.org/10.1016/S0016-7037(03)00273-4.\n#\n")
outFile.write("# Concerning temperature and pressure ranges, the following parameters have been used:\n")
outFile.write(f"# min temperature = {temperatures.min}, max temperature = {temperatures.max}, #temperature sampling points = {temperatures.numSamples}\n")
outFile.write(f"# min phase pressure = {pressures.min}, max phase pressure = {pressures.max}, #pressure sampling points = {pressures.numSamples}\n#\n")
outFile.write("# temperature [°C], phase pressure [Pa],         y_H2O [-],         x_H2 [-]\n")

# In the following, the notation and equation numbers
# are similar to Spycher et al. 2003 and our paper.
for i in range(temperatures.numSamples):
    T = temperatures[i]
    for j in range(pressures.numSamples):
        p = pressures[j]
        rhoH2 = densities[i][j]
        ms = saltMole # salt molality
        # convert temperature to Kelvin for the function calls
        A = computeA(T + 273.15, p, rhoH2); # (11)
        B = computeBsalt(T + 273.15, p, rhoH2, ms); # (12)
        nu = 2
        y_H2O = ((1 - B)*55.508)/((1/A - B)*(nu*ms+55.508)+nu*ms*B) # (13)
        x_H2 = B*(1 - y_H2O) # (14)

        outFile.write(f" {T:.11e},   {p:.11e}, {y_H2O:.11e}, {x_H2:.11e}\n")

print(f"A file {fileName} has been generated.")
