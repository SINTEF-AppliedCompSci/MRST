% UTILS
%
% Files
%   addFluxesFromSourcesAndBC       - Add in fluxes imposed by sources and face boundary conditions
%   assignValue                     - Assign values to ADI object by way of indices, without changing jacobians
%   bc2ADbc                         - INTERNAL DEPRECATED FUNCTION: Intentionally undocumented.
%   checkWellConvergence            - Compute convergence for wells.
%   CNV_MBConvergence               - Compute convergence based on total mass balance and maximum residual mass balance.
%   compressSchedule                - Compress schedule to take the longest possible timesteps while honoring controls
%   convertDeckScheduleToMRST       - Convert deck-type schedule to MRST style schedule
%   convertIncompWellSols           - Convert wellSols from incomp module to format used in ad-core/ad-blackoil
%   convertReportToSchedule         - Create a new schedule based on actual ministeps from a simulation report
%   convertReservoirFluxesToSurface - Compute surface fluxes from reservoir fluxes
%   double2ADI                      - Convert a double to ADI variable, using a sample ADI variable for dimensions
%   fastInterpTable                 - Fast interpolation of table, using griddedInterpolant
%   getBoundaryConditionFluxesAD    - Get boundary condition fluxes for a given set of values
%   getFaceTransmissibility         - Compute face transmissibilities, accounting for input-specific multipliers
%   getMultipliers                  - Get dynamic multiplier values for reservoir quantities
%   getPerforationToWellMapping     - Get map from global perforation number to global well index.
%   getReportMinisteps              - Get the timesteps used for the ministeps of a report
%   getSimulationTime               - Get the global time for a set of states produced by simulateScheduleAD
%   getSourceFluxesAD               - Short description
%   getWellOutput                   - Extract values from wellsols.
%   initWellSolAD                   - Set up well solution struct for a automatic differentiation model
%   makeScheduleConsistent          - Ensure that a schedule is consistent in terms of well counts/perforations
%   padRatesAndCompi                - Pad one/two/threephase values with zeros corresponding to missing phases.
%   pressureBCContrib               - LEGACY FUNCTION: Intentionally undocumented.
%   pressureBCContribADI            - LEGACY FUNCTION: Intentionally undocumented.
%   printConvergenceReport          - Print a neatly formatted convergence report
%   recoverVars                     - Recover previously eliminated variables x at position n using solutions sol
%   refineSchedule                  - Compute a finer schedule, including new time steps but preserving the time steps of the original
%   reorderForILU                   - Attempt to reorder a set of equations so that the diagonal is non-zero
%   ResultHandler                   - Class for storing and retrieving simulation results, either in memory or stored to disk
%   selectModelFromDeck             - Select simulation model from a ECLIPSE/FrontSim style input deck
%   setupOperatorsTPFA              - Set up helper structure for solvers based on automatic differentiation.
%   setWellSign                     - Ensure that wells have a defined sign. Will attempt to guess based on controls.
%   simpleSchedule                  - Make a schedule with varying timesteps and fixed wells/bc/src terms
%   terniaryWellPlot                - Plot well curves (water, gas, oil and optionally BHP) for wellSols

%{
Copyright 2009-2016 SINTEF ICT, Applied Mathematics.

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
