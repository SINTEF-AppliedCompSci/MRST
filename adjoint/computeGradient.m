function grad = computeGradient(W, adjRes, schedule, controls)
% compute gradient for control variables and project according to
% linEqConst A*u = b  => A*grad = 0
% Thus the projected gradient is
%   grad_p  = (I - A'*inv(A*A')*A)*grad

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


bhpWells  = strcmpi('bhp' , { W.type });
rateWells = strcmpi('rate', { W.type });
S         = [ W.S ];
Dw        = blkdiag( S.D );
DwD       = Dw(:, bhpWells);

% collect all linEqConst
ec = controls.linEqConst;
if ~isempty(ec)
    A = ec.A;
    projector = eye( size(A, 2) ) - A'*( (A*A')\A );
end

[A_N, b_N, A_D, b_D] = controls2Wells(W, schedule, controls);

grad = cell([1, numel(A_N)]);
for k = 1 : numel(A_N)
   adjWellPres = [adjRes(k+1).resSol.wellSol.pressure]';
   l_p         = adjWellPres(rateWells);
   l_q         = vertcat(adjRes(k+1).resSol.wellSol.flux);

   grad_k  =  A_N{k}'*l_p + A_D{k}'*DwD'*l_q;   % non-projected gradient
   if ~isempty(ec)
      grad{k}  = projector * grad_k;
   else
      grad{k}  = grad_k;
   end
end

if controls.numControlSteps == 1,
   gradMat = cell2mat(grad);
   avgVal  = mean(gradMat, 2);
   grad(:) = { avgVal };
end
