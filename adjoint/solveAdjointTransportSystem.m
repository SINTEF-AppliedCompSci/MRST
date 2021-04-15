function [adjRes] = solveAdjointTransportSystem(G, S, W, rock, fluid, simRes, adjRes, obj, varargin)

% Find current time step (search for empty slots in adjRes)
% NOTE: actually curent time step +1

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


numSteps = numel(simRes);
if isempty(adjRes)
    curStep = numSteps;
else
    curStep = find( cellfun(@isempty, {adjRes.timeInterval}), 1, 'last');
end
dt      = simRes(curStep).timeInterval * [-1 1]';
adjRes(curStep).timeInterval = simRes(curStep).timeInterval;

% Generate system matrix
numC    = G.cells.num;
PV      = G.cells.volumes.*rock.poro;
invDPV  = spdiags(1./PV, 0, numC, numC);

[mob, dmob] = mobilities(simRes(curStep).resSol, fluid);
Lt          = sum(mob, 2);
f           = bsxfun(@rdivide, mob, Lt);
Dfw_all     = (f(:,2).*dmob(:,1) - f(:,1).*dmob(:,2)) ./ Lt;   % Chain rule.
DDf         = spdiags(Dfw_all, 0, numC, numC);

%DDf     = spdiags( fluid.dfw(simRes(curStep).resSol), 0, numC, numC);

At      = generateUpstreamTransportMatrix(G, S, W, simRes(curStep).resSol, ...
                              simRes(curStep).wellSol, 'Transpose', true);
systMat = speye(numC, numC) - dt * ( DDf * At * invDPV);   % system matrix
clear PV invDPV DDf At

% Generate right-hand-side
RHS     = -obj.partials(curStep).s';
if curStep < numSteps
    numCF   = size(S.B, 1);
    numW    = numel(W);
    RHS     = RHS + adjRes(curStep+1).resSol.s;
    %dim     = fluid.dLtInv(simRes(curStep).resSol);
    dim     = -sum(dmob, 2) ./ Lt.^2;   % d/ds (1/Lt)

    % Update B^{n+1}v^{n+1} part
    cellNo  = rldecode(1:G.cells.num, diff(G.cells.facePos), 2) .';
    sgn     = 2*double(G.faces.neighbors(G.cells.faces(:,1),1) == cellNo) - 1;
    v       = sgn .* simRes(curStep+1).resSol.flux(G.cells.faces(:,1));   % v^{n+1}
    l_v     = sgn .* adjRes(curStep+1).resSol.flux(G.cells.faces(:,1));   % lam_v^{n+1}

    S.C     = sparse(1:numel(cellNo), cellNo, 1);
    RHS     = RHS - S.C'*spdiags( (S.C*dim).*v , 0, numCF, numCF)*S.B*l_v;

    % Update B_w^{n+1}q_w^{n+1} part
    for wellNr = 1:numW
        w   = W(wellNr);
        q   = simRes(curStep+1).resSol.wellSol(wellNr).flux; % q_w^{n+1}
        l_q = adjRes(curStep+1).resSol.wellSol(wellNr).flux; % lam_q_w^{n+1}
        RHS = RHS + w.S.C'*( ( (w.S.C*dim).*q ).*(w.S.B*l_q) );
    end
end

% Solve system
adjRes(curStep).resSol.s  = systMat \ RHS;

end

%--------------------------------------------------------------------------

function [mob, dmob] = mobilities(state, fluid)
   %output/derivatives should be wrt s_w
   mu = fluid.properties(state);
   s  = fluid.saturation(state);
   [kr{1:2}] = fluid.relperm(s, state);

   %        \lambda_i in varargout{1}.
   % (d/ds) \lambda_i in varargout{2}.  Returned iff requested.
   %
   mob = bsxfun(@rdivide, kr{1}, mu);
   if nargout > 1
       dmob = bsxfun(@rdivide, kr{2}(:, [1 end]), mu);
       dmob(:, 2) = -dmob(:,2);
   end
   %kr = cellfun(@(x)x(:,[1 end]), kr, 'UniformOutput', false);
   %varargout = cellfun(@(n) bsxfun(@rdivide, n, mu), kr, ...
   %                    'UniformOutput', false);
end

%--------------------------------------------------------------------------

% function varargout = mobilities(state, fluid)
%    mu = fluid.properties(state);
%    s  = fluid.saturation(state);
%    [kr{1:nargout}] = fluid.relperm(s, state);
%
%    %        \lambda_i in varargout{1}.
%    % (d/ds) \lambda_i in varargout{2}.  Returned iff requested.
%    %
%    kr = cellfun(@(x)x(:,[1 end]), kr, 'UniformOutput', false);
%    varargout = cellfun(@(n) bsxfun(@rdivide, n, mu), kr, ...
%                        'UniformOutput', false);
% end
