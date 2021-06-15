% UTILS
%
% Files
%   addPropertyDependence                 - Document dependencies and external dependencies
%   addFluxesFromSourcesAndBC             - Add in fluxes imposed by sources and face boundary conditions
%   assignValue                           - Assign values to ADI object by way of indices, without changing jacobians
%   assignWellValuesFromControl           - Assign wellSol values when values are set as controls
%   bc2ADbc                               - INTERNAL DEPRECATED FUNCTION: Intentionally undocumented.
%   checkWellConvergence                  - Compute convergence for wells.
%   CNV_MBConvergence                     - Compute convergence based on total mass balance and maximum residual mass balance.
%   calculatePhaseRate                    - Undocumented Utility Function
%   combineEquations                      - Combine equations. For doubles, this is equivialent to a vertical
%   combineSchedules                      - Combine multiple schedules to form a schedule with multiple controls
%   compressSchedule                      - Compress schedule to take the longest possible timesteps while honoring controls
%   computeCpGeometry                     - Undocumented Utility Function
%   computeSourcesAndBoundaryConditionsAD - Compute phase-pseudocomponent source terms (compatible with AD codes)
%   convert2MSWell                        - Utility for Converting Standard Well Structure to Multi-Segment Type
%   convertDeckScheduleToMRST             - Convert deck-type schedule to MRST style schedule
%   convertIncompWellSols                 - Convert wellSols from incomp module to format used in ad-core/ad-blackoil
%   convertReportToSchedule               - Create a new schedule based on actual ministeps from a simulation report
%   convertReservoirFluxesToSurface       - Compute surface fluxes from reservoir fluxes
%   criticalPointChop                     - Perform one or two-sided stability chop for an updated value
%   crossFlowMixture                      - Undocumented Utility Function
%   crossFlowMixtureDensity               - Undocumented Utility Function
%   double2ADI                            - Convert a double to ADI variable, using a sample ADI variable for dimensions
%   estimateCompositionCFL                - Undocumented Utility Function
%   estimateSaturationCFL                 - Undocumented Utility Function
%   expandIfUniform                       - Utility which reverses "value" compaction. If given a matrix (logical
%   expandMatrixToCell                    - Expand a matrix into cell arrays. Typical usage: Converting state
%   faceUpstr                             - Perform single-point upwinding of cell values to face
%   fastInterpTable                       - Fast interpolation of table, using griddedInterpolant
%   filterSchedule                        - Filter unused controls from a schedule
%   getBoundaryConditionFluxesAD          - Get boundary condition fluxes for a given set of values
%   getCellMajorReordering                - Get equation ordering transforming variable major to cell major ordering
%   getConvergenceValuesCNV               - Compute convergence based on total mass balance and maximum residual mass balance.
%   getConvergenceValuesWells             - Undocumented Utility Function
%   getEquilPC                            - Undocumented Utility Function
%   getFractionalFlowMagnitude            - Undocumented Utility Function
%   getGridSYMRCMOrdering                 - Undocumented Utility Function
%   getMultiDimInterpolator               - Get a multidimensional interpolator (with support for ADI varibles)
%   getMultipliers                        - Get dynamic multiplier values for reservoir quantities
%   getPerforationToWellMapping           - Get map from global perforation number to global well index.
%   getReportMinisteps                    - Get the timesteps used for the ministeps of a report
%   getReportOutput                       - Get output from report after call to simulateScheduleAD
%   getReservoirModel                     - Get the underlying reservoir model of a WrapperModel
%   getSampleAD                           - Utility for getting a AD value if it exists from a list of possible
%   getSimulationTime                     - Get the global time for a set of states produced by simulateScheduleAD
%   getSourceFluxesAD                     - Short description
%   getWellOutput                         - Extract values from wellsols.
%   HandleStruct                          - 
%   initWellSolAD                         - Set up well solution struct for a automatic differentiation model
%   interpolateIDW                        - Undocumented Utility Function
%   makeScheduleConsistent                - Ensure that a schedule is consistent in terms of well counts/perforations
%   mergeOrderedArrays                    - Merge two sets of cells that are similar in that they may contain
%   numelData                             - Alias for numel. Useful for writing code which handles either
%   numelValue                            - Undocumented Utility Function
%   padRatesAndCompi                      - Pad one/two/threephase values with zeros corresponding to missing phases.
%   phaseDensitiesTobfactor               - Convert densities to b-facctors, accounting for dissolution
%   pressureBCContrib                     - LEGACY FUNCTION: Intentionally undocumented.
%   pressureBCContribADI                  - LEGACY FUNCTION: Intentionally undocumented.
%   printConvergenceReport                - Print a neatly formatted convergence report
%   readSummaryLocal                      - Undocumented Utility Function
%   recoverVars                           - Recover previously eliminated variables x at position n using solutions sol
%   refineSchedule                        - Compute a finer schedule, including new time steps but preserving the time steps of the original
%   reorderForILU                         - Attempt to reorder a set of equations so that the diagonal is non-zero
%   ResultHandler                         - Class for storing and retrieving simulation results, either in memory or stored to disk
%   selectLinearSolverAD                  - Undocumented Utility Function
%   selectModelFromDeck                   - Select simulation model from a ECLIPSE/FrontSim style input deck
%   setupOperatorsTPFA                    - Set up helper structure for solvers based on automatic differentiation.
%   setMPFADiscretization                 - Set MPFA discretization on a model
%   setReservoirModel                     - Set the underlying reservoir model of a WrapperModel
%   setTimeDiscretization                 - Set the discretization choice for a model
%   setWellSign                           - Ensure that wells have a defined sign. Will attempt to guess based on controls.
%   setWENODiscretization                 - Set WENO discretization on a model
%   simpleSchedule                        - Make a schedule with varying timesteps and fixed wells/bc/src terms
%   sizeJac                               - Undocumented Utility Function
%   splitFaceCellValue                    - Split multi-valued function into cell and face values
%   splitMatrixForReduction               - Split matrix A and right-hand side into blocks
%   standaloneSolveAD                     - Solve a single time-step with AD solvers for given forces
%   structPropEvaluated                   - Undocumented Utility Function
%   terniaryWellPlot                      - Plot well curves (water, gas, oil and optionally BHP) for wellSols
%   wellSolToVector                       - Extract selected summary vectors from cell array of well solutions

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
