% Files
%   SolveEqsADI.m              - We have been provided a pod basis, use it to reduce equation sets
%   computeNumGrad.m           - Compute numerical gradient w.r.t wells
%   dinterpReg.m               - if reginx{k} == ':' for improved eff seperate this case
%   getConvergence.m           - Compute convergence based on total mass balance (tol_mb) and maximum
%   getRegMap.m                - if nargin < 4 whole domain
%   getResiduals.m             - Store the residuals for debugging and convergence testing.
%   getSimMatrices.m           - half-trans -> trans and reduce to interior
%   getWellStuff.m             - ------------------------------------------
%   handleRegions.m            - also need actnum/G.cells.indexmap
%   initWellSolLocal.m         - Initialize well solution data structure.
%   interpReg.m                - if reginx{k} == ':' for improved eff seperate this case
%   printResidual.m            - fprintf('-9s', eqnnames{:})
%   scheduleFromSummary.m      - There are n-1 control steps for n summary steps
%   setupSimComp.m             - Set up helper structure for solvers based on automatic differentiation.
%   solveAdjointEqsADI.m       - If adjVec is not empty, we are not at the last timestep (first in the
%   updateSchedule.m           - --------------------------------------------------------------------------
%   wellSolToVector.m          - Helper function which makes cell arrays of well solutions easier to plot

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
