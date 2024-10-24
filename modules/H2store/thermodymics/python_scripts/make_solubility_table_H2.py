#!/usr/bin/env python3

""" Generate table for H2 - H2O solubilities.

The calculations implement model from
Raad et al. (2022),
https://doi.org/10.1016/j.ijhydene.2022.09.222.
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


def equilibriumConstantH2(T):
    TinC = T - 273.15 # temperature in °C
    c = [2.9947, 4.81373e-3, -5.1773e-5, 1.19052e-7]
    logk0_H2 = c[0] + c[1] * TinC + c[2] * TinC * TinC + c[3] * TinC * TinC*TinC
    return math.pow(10, logk0_H2)

def equilibriumConstantH2O(T):
    TinC = T - 273.15 # temperature in °C
    c = [-2.1817, 2.98e-2, -1.098e-4, 2.048e-7]
    logk0_H2O = c[0] + c[1]*TinC + c[2]*TinC*TinC + c[3]*TinC*TinC*TinC
    return math.pow(10, logk0_H2O)   

def fugacityCoefficientH2(T, p, rhoH2):
    molarMassH2 = 2.016e-3 # [kg/mol]
    V = 1/(rhoH2/molarMassH2)*1e6 # molar volume [cm3/mol]
    p_bar = p/1e5 # phase pressure in bar
    a_H2 = 1441753.379 # mixture parameter of Redlich-Kwong equation
    b_H2 = 18.417 # mixture parameter of Redlich-Kwong equation
    R = 83.1446261815324 # universal gas constant [bar.cm3/mol.K]

    lnPhiH2 = math.log(V/(V - b_H2)) + b_H2/(V - b_H2) \
               - 2*a_H2/(R*math.pow(T, 1.5)*b_H2)*math.log((V + b_H2)/V) \
               + a_H2*b_H2/(R*math.pow(T, 1.5)*b_H2*b_H2) \
                 *(math.log((V + b_H2)/V) - b_H2/(V + b_H2)) \
               - math.log(p_bar*V/(R*T))
    return math.exp(lnPhiH2)

def fugacityCoefficientH2O(T, p, rhoH2):
    molarMassH2 = 2.016e-3 # [kg/mol]
    V = 1/(rhoH2/molarMassH2)*1e6 # molar volume [cm3/mol]
    p_bar = p/1e5 # phase pressure in bar
    a_H2 = 1441753.379 # mixture parameter of  Redlich-Kwong equation
    a_H2O = 142666655.8 # mixture parameter of  Redlich-Kwong equation
    a_H2_H2O = math.pow(a_H2*a_H2O,0.5) # mixture parameter of Redlich-Kwong equation
    b_H2 = 18.417 # mixture parameter of Redlich-Kwong equation
    b_H2O = 21.127 # mixture parameter of Redlich-Kwong equation
    R = 83.1446261815324 # universal gas constant [bar.cm3/mol.K]

    lnPhiH2O = math.log(V/(V - b_H2)) + b_H2O/(V - b_H2) \
               - 2*a_H2_H2O/(R*math.pow(T, 1.5)*b_H2)*math.log((V + b_H2)/V) \
               + a_H2*b_H2O/(R*math.pow(T, 1.5)*b_H2*b_H2) \
                 *(math.log((V + b_H2)/V) - b_H2/(V + b_H2)) \
               - math.log(p_bar*V/(R*T))
    return math.exp(lnPhiH2O)

def computeA(T, p, rhoH2):
    deltaP = p/1e5 - 1 # pressure range [bar] from p0 = 1 bar to p
    v_av_H2O = 18.1 # average partial molar volume of H2O [cm3/mol]
    R = 83.1446261815324 # universal gas constant [bar.cm3/mol.K]
    k0_H2O = equilibriumConstantH2O(T) # equilibrium constant for H2O at 1 bar
    phi_H2O = fugacityCoefficientH2O(T, p, rhoH2) # fugacity coefficient of H2O for the water-H2 system
    p_bar = p/1e5
    return k0_H2O/(phi_H2O*p_bar)*math.exp(deltaP*v_av_H2O/(R*T))

def computeB(T, p, rhoH2):
    deltaP = p/1e5 - 1 # pressure range [bar] from p0 = 1 bar to p
    v_av_H2 = 26.7 # average partial molar volume of H2 [cm3/mol]
    R = 83.1446261815324 # universal gas constant [bar.cm3/mol.K]
    k0_H2 = equilibriumConstantH2(T) # equilibrium constant for H2 at 1 bar
    phi_H2 = fugacityCoefficientH2(T, p, rhoH2) # fugacity coefficient of H2 for the water-H2 system
    p_bar = p/1e5
    return phi_H2*p_bar/(55.508*k0_H2)*math.exp(-(deltaP*v_av_H2)/(R*T))


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
cmdArgs = vars(parser.parse_args())

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
outFile = open(fileName, "w")
outFile.write(f"# This autogenerated file contains solubilities of H2 and H2O in a respective fluid system.\n")
outFile.write("# The values have been calculated by means of (11)-(14) in https://doi.org/10.1016/S0016-7037(03)00273-4.\n#\n")
outFile.write("# Concerning temperature and pressure ranges, the following parameters have been used:\n")
outFile.write(f"# min temperature = {temperatures.min}, max temperature = {temperatures.max}, #temperature sampling points = {temperatures.numSamples}\n")
outFile.write(f"# min phase pressure = {pressures.min}, max phase pressure = {pressures.max}, #pressure sampling points = {pressures.numSamples}\n#\n")
outFile.write("# temperature [°C], phase pressure [Pa],         y_H2O [-],         x_H2 [-]\n")

# In the following, the notation and equation numbers
# from Spycher et al. 2003 are used.
for i in range(temperatures.numSamples):
    T = temperatures[i]
    for j in range(pressures.numSamples):
        p = pressures[j]
        rhoH2 = densities[i][j]
        # convert temperature to Kelvin for the function calls
        A = computeA(T + 273.15, p, rhoH2); # (11)
        B = computeB(T + 273.15, p, rhoH2); # (12)
        y_H2O = (1 - B)/(1/A - B) # (13)
        x_H2 = B*(1 - y_H2O) # (14)

        outFile.write(f" {T:.11e},   {p:.11e}, {y_H2O:.11e}, {x_H2:.11e}\n")

print(f"A file {fileName} has been generated.")
