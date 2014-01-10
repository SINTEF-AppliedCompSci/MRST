% ADJOINT Add-on module for adjoint-based production optimisation
%
% Files
%   addAdjointWellFields              - SYNOPSIS:
%   adjointFluidFields                - Extend fluid functionality with fields needed in (2nd order) adjoint imp.
%   assembleWellSystem                - Generate pressure linear system components for wells.
%   computeAdjointRHS                 - Compute adjoint 'pressure' rhs
%   computeGradient                   - compute gradient for control variables and project according to
%   computeNumericalGradient          - compute numerical gradient
%   computeNumericalGradientMS        - compute numerical gradient
%   controls2RHS                      - Create mappings A_N, b_N, A_D, b_D such that
%   controls2Wells                    - Create mappings A_N, b_N, A_D, b_D such that
%   dispControls                      -
%   dispSchedule                      -
%   generateNonUniformCoarseGridv2    - Hacked from Veras version .....
%   generateSat2MS                    - Create mappings from coarse saturation field to MS - pressure system
%   generateUpstreamTransportMatrix   - generateUpstreamTransportMatrix for use in saturation solver
%   generateUpstreamTransportMatrixMS - generateUpstreamTransportMatrixMS for use in saturation solver
%   initControls                      - initControls -- Initialize control structure based on well schedule
%   initSchedule                      - initSchedule -- Initialize schedule structure based on well W.
%   lineSearch                        - Run ad-hoc line search based on given gradient.
%   lineSearch2                       - Run ad-hoc line search based on given gradient
%   lineSearchMS                      - Run ad-hoc line search based on given gradient
%   optimizeObjective                 - optimizeObjective -- Run whole optimization proccess
%   optimizeObjectiveMS               - optimizeObjectiveMS -- Run whole optimization proccess multiscale
%   projectGradient                   - Project gradient according to linear input constraints. Handles box-constraints and
%   runAdjoint                        - runAdjoint -- Run adjoint simulation based on simRes and schedule.
%   runAdjointMS                      - runAdjoint -- Run adjoint simulation based on simRes and schedule.
%   runSchedule                       - runSchedule -- Run simulation based on schedule.
%   runScheduleMS                     - runSchedule -- Run simulation based on schedule.
%   solveAdjointPressureSystem        - Find current time step (search for empty slots in adjRes)
%   solveAdjointPressureSystemMS      - Find current time step (search for empty slots in adjRes)
%   solveAdjointTransportSystem       - Find current time step (search for empty slots in adjRes)
%   solveAdjointTransportSystemMS     - DLtInv = @(sol)(-fluid.dkr(sol)*(1./fluid.mu)'...
%   updateSchedule                    - Update schedule based on controls
%   updateWells                       - Update wells based on schedule time step
