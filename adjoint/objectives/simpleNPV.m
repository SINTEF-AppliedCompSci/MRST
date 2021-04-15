function obj = simpleNPV(G, S, W, rock, fluid, simRes, ...
                         schedule, controls, varargin)
%Simple net-present-value function - no discount factor
%
% SYNOPSIS:
%   obj = simpleNPV(G, S, W, rock, fluid, simRes, schedule, controls)
%
% DESCRIPTION:
%   Computes value of objective function for given simulation, and partial
%   derivatives of variables if varargin > 6
%
% PARAMETERS:
%   simRes      -
%
% RETURNS:
%   obj         - structure with fields
%
% SEE ALSO:

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


opt     = struct('OilPrice'              , 100   , ...
                 'WaterProductionCost'   ,  10 , ...
                 'WaterInjectionCost'    ,  10 , ...
                 'DiscountFactor', 0     , ...
                 'ComputePartials',         []);
opt     = merge_options(opt, varargin{:});

ro      = opt.OilPrice            / stb;
rw      = opt.WaterProductionCost / stb;
ri      = opt.WaterInjectionCost  / stb;
d       = opt.DiscountFactor;

%-----------------------------------------------
computePartials = opt.ComputePartials;
if isempty(computePartials)
    computePartials  = (nargin > 6);
end
numSteps = numel(simRes);
val      = 0;
partials = repmat( struct('v', [], 'p', [], 'pi', [], ...
                          's', [], 'u', []), [numSteps 1] );

totTime  = max( [simRes.timeInterval] );

mob = cell([1, 1 + 2*double(computePartials)]);

for step = 2 : numSteps,
    resSol  = simRes(step).resSol;
    wellSol = simRes(step).wellSol;
    int     = simRes(step).timeInterval;
    dt      = int(2) - int(1);
    dFac    = (1+d)^(-int(2)/year);

    [wellRates, rateSigns] = getRates(W, wellSol);
    wellCells = vertcat( W.cells );
    wellSats  = resSol.s( wellCells );

    [mob{:}] = mobilities(struct('s', wellSats), fluid);

    Lt  = sum(mob{1}, 2);
    f   = bsxfun(@rdivide, mob{1}, Lt);

    f_w = f(:,1);
    f_o = 1 - f_w;

    injInx  = (rateSigns > 0);
    prodInx = (rateSigns < 0);

    % Objective value:

    val   = val + dt*dFac*( - sum(  wellRates(injInx)                )*ri ...
                            - sum( -wellRates(prodInx).*f_w(prodInx) )*rw ...
                            + sum( -wellRates(prodInx).*f_o(prodInx) )*ro );

    if computePartials,
        numCF    = size(G.cells.faces, 1);
        numC     = G.cells.num;
        numF     = G.faces.num;
        numU     = numel(controls.well);

        partials(step).v   = zeros(1, numCF);
        partials(step).p   = zeros(1, numC);
        partials(step).pi  = zeros(1, numF);

        partials(step).q_w =  - dt*dFac*ri*injInx' ...
                              + dt*dFac*rw*( prodInx.*f_w )' ...
                              - dt*dFac*ro*( prodInx.*f_o )';

        df  = (f(:,2).*mob{2}(:,1) - f(:,1).*mob{2}(:,2)) ./ Lt;   % Chain rule.
        d2f = (f(:,2).*mob{3}(:,1) - f(:,1).*mob{3}(:,2) - ...
               2 .* df .* sum(mob{2}, 2)) ./ Lt;

        Df_w = df(:,1);
        Df_o = - Df_w;

        ds   = zeros(1, numC);
        ds( wellCells(prodInx) )  =  dt*dFac*rw*( wellRates(prodInx).*Df_w(prodInx) )' ...
                                    -dt*dFac*ro*( wellRates(prodInx).*Df_o(prodInx) )';

        partials(step).s = ds;
        partials(step).u = zeros(1, numU);

        % Disable this...
        %%{
        % Second order derivatives:
        %if isfield(fluid, 'd2fw')
        %D2f_w = fluid.d2fw( struct('s', wellSats) );
        D2f_w = d2f(:,1);
        D2f_o = -D2f_w;

        d2s = zeros(numC, 1);
        d2s( wellCells(prodInx) )  =  dt*dFac*rw*( wellRates(prodInx).*D2f_w(prodInx) )' ...
           -dt*dFac*ro*( wellRates(prodInx).*D2f_o(prodInx) )';
        partials(step).s2   = spdiags(d2s, 0, numC, numC);

        % --------------------
        %IX = sparse( wellCells(prodInx), (1:nnz(prodInx))', 1, G.cells.num, nnz(prodInx));
        % ds' = IX * (dt*dFac*rw*prodRates*(IX' * Df_w(s)) + ...
        %dssl = dt*dFac*rw*( wellRates(prodInx).*D2f_w(prodInx) ) ...
        %      -dt*dFac*ro*( wellRates(prodInx).*D2f_o(prodInx) );
        %partss = IX*spdiags(dssl, 0,  nnz(prodInx), nnz(prodInx))*IX';

        partials(step).qs   = dt * dFac * sparse(wellCells(prodInx), find(prodInx), rw*Df_w(prodInx)-ro*Df_o(prodInx), ...
           numC, length(prodInx));
        %end
        %}
    end
end

obj.val = val;
if computePartials, obj.partials = partials; end
end

%--------------------------------------------------------------------------

function [wellRates, rateSigns] = getRates(W, wellSol)
   wellRates = vertcat(wellSol.flux);
   wellSys   = [ W.S ];
   Dw        = blkdiag( wellSys.D );
   wellSigns = ones( numel(W), 1 );
   totRates  = Dw'*wellRates;
   wellSigns( totRates < 0 ) = -1;

   for k = 1:numel(W),
      if ~isempty(W(k).sign) % override if sign is given expl
         wellSigns(k) = W(k).sign;
      end
   end

   rateSigns = Dw*wellSigns;
end

%--------------------------------------------------------------------------

function [mob, dmob, dmob2] = mobilities(state, fluid)
   %output/derivatives should be wrt s_w
   mu = fluid.properties(state);
   s  = fluid.saturation(state);
   [kr{1:nargout}] = fluid.relperm(s, state);

   %        \lambda_i in varargout{1}.
   % (d/ds) \lambda_i in varargout{2}.  Returned iff requested.
   %
   mob = bsxfun(@rdivide, kr{1}, mu);
   if nargout > 1
       dmob = bsxfun(@rdivide, kr{2}(:, [1 end]), mu);
       dmob(:, 2) = -dmob(:,2);
   end
    if nargout > 2
       dmob2 = bsxfun(@rdivide, kr{3}(:, [1 end]), mu);
       %dmob2(:, 2) = -dmob(:,2);
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

%{
function [wellRates, rateSigns] = getRates(W, wellSol)
wellRates   = vertcat(wellSol.flux);
wellSys     = [ W.S ];
Dw    = blkdiag( wellSys.D );
if isfield(W, 'sign')
    wellSigns = vertcat(W.sign);
else
    wellSigns    = ones( numel(W), 1 );
    totRates = Dw'*wellRates;
    wellSigns( totRates < 0 ) = -1;
end
%}
