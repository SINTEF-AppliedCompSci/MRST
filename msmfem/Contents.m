% Routines supporting the multiscale mixed FE method for the pressure equation.
%
% Files
%   basisMatrixHybrid            - Form hybrid versions of the multiscale basis function matrices.
%   basisMatrixMixed             - Form mixed versions of the multiscale basis function matrices.
%   dynamicCoarseWeight          - Compute synthetic multiscale weighting function.
%   generateCoarseSystem         - Construct coarse system component matrices from fine grid model.
%   generateCoarseWellSystem     - Construct coarse system component matrices for well contributions.
%   solveCoarsePsysBO            - Solve coarsened fine-scale well system (for Black Oil).
%   solveIncompFlowMS            - Solve coarse (multiscale) pressure system.
%   solveIncompFlowMSSpeedUp     - Solve (multiscale) pressure system assuming that some values are precomputed.
%   speedUpMS                    - Precompute MsMFE basis reduction matrices in order speed up assembly
%   unpackWellSystemComponentsMS - Extract coarse linear system components from wells.
%   updateBasisFunc              - Update basis functions in regions where the total mobility has changed

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
