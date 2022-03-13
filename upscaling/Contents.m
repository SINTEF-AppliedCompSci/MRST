% UPSCALING Solvers for permeability upscaling problem
%
% Files
%   computeMimeticIPGp          - Wrapper for computeMimeticIP for periodic grids.
%   computeTransGp              - Wrapper for computeTrans for periodic grids.
%   makeMatchingGridFromGrdel   - {
%   makePeriodicGridMulti3d     - Construct a periodic grid.
%   simulateToSteadyState       - Simulate until steady state is reached
%   upscaleGridRockAndWells     - [rock_cg,CG,W_cg] = upscaleGridRockAndWells(G,rock,coarseDim)
%   upscalePC                   - {
%   upscalePCCaplimit           - {
%   upscalePerm                 - Compute upscaled permeabilites using flow based upscaling.
%   upscalePermeabilityFixed    - function for doing permeability upscaling on grids with fixed pressure conditions (pressure drop in one direction and noflow elsewhere)
%   upscalePermeabilityMim      - function for doing permeability upscaling on periodic grids
%   upscalePermeabilityPeriodic - function for doing permeability upscaling on periodic grids
%   upscaleRelperm              - Upscale relative permeability.
%   upscaleRelpermLimit         - Upscale relperm based on viscous/capillary limit tabulated by saturation
%   upscaleTrans                - Calculate upscaled transmissibilities for a coarse model

%{
Copyright 2009-2022 SINTEF Digital, Mathematics & Cybernetics.

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
