function [adjRes] = solveAdjointPressureSystem(G, S, W, rock, fluid, simRes, adjRes, obj, varargin)

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


curStep = find( cellfun(@(x)~isempty(x), {adjRes.timeInterval}), 1, 'first');
dt      = simRes(curStep).timeInterval * [-1 1]';

% Generate RHS, that is f-part (rest is zero)
PV      = G.cells.volumes.*rock.poro;
invPV   = 1./PV;
%f_w     = fluid.krw(simRes(curStep).resSol) ./
%fluid.Lt(simRes(curStep).resSol); % fractinal flow, f(s^n)
%f_w     = fluid.fw( simRes(curStep).resSol);

%{
mu  = fluid.properties(simRes(curStep).resSol);
sat = fluid.saturation(simRes(curStep).resSol);
kr  = fluid.relperm(sat, simRes(curStep).resSol);
mob = bsxfun(@rdivide, kr, mu);
f_w = mob(:,1) ./ sum(mob,2);
%}
[mob, dmob] = mobilities(simRes(curStep).resSol, fluid);
Lt          = sum(mob, 2);
f           = bsxfun(@rdivide, mob, Lt);
f_w         = f(:,1);

l_s     = adjRes(curStep).resSol.s;        % lam_s^n
% Flux-matrix: A.i, A.j, A.qMinus
[A, qPluss, signQ] = generateUpstreamTransportMatrix(G, S, W, simRes(curStep).resSol, ...
                                        simRes(curStep).wellSol, 'VectorOutput', true);

dQPluss =  double( signQ > 0 );
dQMinus = -double( signQ < 0 );

f_bc    = -obj.partials(curStep).v';

cellNo  = rldecode(1:G.cells.num, diff(G.cells.facePos), 2) .';
S.C     = sparse(1:numel(cellNo), cellNo, 1);
f_bc    = f_bc - dt*( f_w(A.j).*( S.C*(invPV.*l_s) ) ) ...
               + dt*( S.C*( (-f_w.*dQMinus + dQPluss).*( invPV.*l_s ) ) );
%f_bc    = f_bc - dt*( f_w(A.j).*( S.C*(invPV.*l_s) ) ) ...
%               + dt*( S.C*( (f_w.*s_qm + s_qp).*( invPV.*l_s ) ) );

% Set f and h in W(i).S.RHS equal DJ/Dp_w and zero. Leave naumannFaces and dirichletFaces as is
inx = 0;
for wellNr = 1 : numel(W)
    numCells = length( W(wellNr).cells );
    W(wellNr).S.RHS.f = obj.partials(curStep).q_w( inx+1 : inx+numCells )';
    W(wellNr).S.RHS.h = zeros( size(W(wellNr).S.RHS.h) );
    inx = inx + numCells;
end

% Solve linear system based on s^{n-1}
b = computeAdjointRHS(G, W, f_bc);

if strcmp(S.type, 'hybrid')
   solver = 'hybrid';
else
   solver = 'mixed';
end

% Solve linear system based on s^{n-1}
resSol = solveIncompFlowLocal(simRes(curStep-1).resSol, G, S, fluid, ...
                         'wells', W, 'rhs', b, 'Solver', solver);

% Update adjRes !!! Note minuses in front of pressure and wellrates in
% forward system, but not in adjoint, thus set minus here.
% REMEMBER: adjRes(curStep).resSol.s has been solved already by
% solveAdjointTransportSystem. Therefore, do not set
% adjRes(curStep).resSol = resSol;
% because this will overwrite the correct value.

adjRes(curStep).resSol.flux = resSol.flux;
adjRes(curStep).resSol.pressure = - resSol.pressure; % Negative!!
adjRes(curStep).resSol.wellSol = resSol.wellSol;
for k = 1 : numel(adjRes(curStep).resSol.wellSol),
    adjRes(curStep).resSol.wellSol(k).flux = ...
       - adjRes(curStep).resSol.wellSol(k).flux;     % !!! Negative
end
end
%{
adjRes(curStep).resSol.cellPressure = - resSol.pressure;           % !!!minus
adjRes(curStep).resSol.facePressure = resSol.facePressure;
adjRes(curStep).wellSol             = resSol.wellSol;
%}

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

% function varargout = mobilities(state, fluid)
%    mu = fluid.properties(state);
%    s  = fluid.saturation(state);
%    [kr{1:nargout}] = fluid.relperm(s, state);
%
%    %        \lambda_i in varargout{1}.
%    % (d/ds) \lambda_i in varargout{2}.  Returned iff requested.
%    %
%     kr = cellfun(@(x)x(:,[1 end]), kr, 'UniformOutput', false);
%    varargout = cellfun(@(n) bsxfun(@rdivide, n, mu), kr, ...
%                        'UniformOutput', false);
% end
