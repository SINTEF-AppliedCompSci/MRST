% UTILS
%
% Files
%   checkComponentMassBalance         - Check mass balance of a simulator run and print to screen
%   coolPropFluidsStructs             - Fluids taken from CoolProp. Used for lookup of properties.
%   cubicPositive                     - Straightforward implementation of a cubic root solver for vectorized
%   ensureMinimumFraction             - Set a minimum value on a composition matrix
%   equationsCompositional            - Overall composition fully-implicit equations
%   estimateEquilibriumWilson         - Estimate equilibrium constant for a given pressure and temperature
%   expandMatrixToCell                - Expand a matrix into cell arrays. Typical usage: Converting state
%   FastAD                            - A very limited AD class for quick EOS assembly
%   formatMassString                  - Small utility which returns a human readable string from mass.
%   getComponentsTwoPhaseSimpleWater  - Undocumented Utility Function
%   getDefaultFlashNonLinearSolver    - Get default nonlinear solver for flash problems.
%   getEOSComponent                   - Undocumented Utility Function
%   getImpesWeightsOverallComposition - Undocumented Utility Function
%   getNonUnitMassFraction            - Internal utility. Intentionally undocumented.
%   getNonUnitMoleFraction            - Internal utility. Intentionally undocumented.
%   getPartialVolumes                 - Undocumented Utility Function
%   initCompositionalState            - Initialize a compositional state given initial composition
%   initDeckEOSModel                  - Set up a EquationOfState model from a parsed deck
%   initVariablesFastAD               - Initialise FastAD objects used internally in some compositional code
%   mrstCubic                         - Straightforward implementation of a cubic root solver for vectorized
%   newtonFugacityEquilibrium         - Single step of the newton update for flash equations
%   phaseStabilityTest                - Perform a phase stability test for a mixture
%   solveRachfordRiceVLE              - Solve Rachford Rice equations to find liquid and vapor
%   standaloneFlash                   - Utility for flashing without explicitly forming a state
%   validateCompositionalForces       - Undocumented Utility Function

%{
Copyright 2009-2021 SINTEF Digital, Mathematics & Cybernetics.

This file is part of The MATLAB Reservoir Simulation Toolbox (MRST).

MRST is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

MRST is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with MRST.  If not, see <http://www.gnu.org/licenses/>.
%}
