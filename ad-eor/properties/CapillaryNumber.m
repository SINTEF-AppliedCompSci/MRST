classdef CapillaryNumber < StateFunction
    
    properties
        vmeth
    end
    
    methods
        function prop = CapillaryNumber(model, varargin)
            prop@StateFunction(model, varargin{:});
            prop = prop.dependsOn({'pressure', 'surfactant', 'wellSol'}, 'state');
            assert(model.water && isfield(model.fluid,'ift'));
            % We can choose square or linear velocity method here.
            prop.vmeth = 'square';
        end
        
        function Nc = evaluateOnDomain(prop, model, state)
            velocitymethod = prop.vmeth;
            add_well_contrib = ~isempty(model.FacilityModel.WellModels);
            if add_well_contrib
                wellSol = model.getProps(state, 'wellSol');
            else
                wellSol = [];
            end
            [p, cs] = model.getProps(state, 'pressure', 'surfactant');
            fluid = model.fluid;
            G = model.G;
            s = model.operators;
            if add_well_contrib
                drivingForces = model.getValidDrivingForces();
                W = drivingForces.W;
                [wellVars, wellVarNames, wellMap] = ...
                    model.FacilityModel.getAllPrimaryVariables(wellSol);
                pBH = wellVars{wellMap.isBHP};
            end
            
            % If the injection well pressure or flow rate is set too large,
            % the pressure change may reach the upper limit at the first time step and
            % cause the pressure gradient (gradp) to be 0 at next simulation step.
            % This may further cause the calculation error of log(Nc) in
            % SurfactantRelativePermeability function. Here we set gradp to be a
            % minimum (1e-8) to avoid such numerical calculation errors.
            gradp = s.Grad(p);
            tooSmall = abs(value(gradp)) < 1e-8;
            gradp = ~tooSmall.*gradp + tooSmall.*1e-8;
            v = -s.T.*gradp;
            
            switch velocitymethod
                case 'linear'
                    veloc = s.veloc;
                    veloc_sq = 0;
                    for i = 1 : numel(veloc)
                        veloc_sq = veloc_sq + veloc{i}(v).^2;
                    end
                    if add_well_contrib
                        [velocW, wc] = computeWellContrib(G, W, p, pBH);
                        veloc_sq(wc) = veloc_sq(wc) + velocW.^2;
                    end                    
                case 'square'
                    veloc_sq = s.sqVeloc(v);
                    % add_well_contrib = false;
                    % if add_well_contrib
                    %    [velocW, wc] = computeWellContrib(G, W, p, pBH);
                    %    veloc_sq(wc) = veloc_sq(wc) + velocW.^2;
                    % end                  
                otherwise
                    error('option for velocCompMethod not recognized');
            end            
            abs_veloc = (veloc_sq).^(1/2);
            sigma     = model.fluid.ift(cs);
            Nc        = abs_veloc./sigma;
        end        
    end
end

% ------------------------------------------------------------------------
function [velocW, wc] = computeWellContrib(G, W, p, pBH)
% We use the bottom hole pressure to compute well velocity
perf2well = getPerforationToWellMapping(W);
Rw = sparse(perf2well, (1:numel(perf2well))', 1, numel(W), numel(perf2well));

wc  = vertcat(W.cells);
pW  = p(wc);
pBH = Rw'*pBH;
Tw  = vertcat(W.WI);
Tw  = Rw'*Tw;
welldir = { W.dir };
i = cellfun(@(x)(numel(x)), welldir) == 1;
welldir(i) = arrayfun(@(w) repmat(w.dir, [ numel(w.cells), 1 ]), ...
    W(i), 'UniformOutput', false);
welldir = vertcat(welldir{:});
[dx, dy, dz] = cellDims(G, wc);
thicknessWell = dz;
thicknessWell(welldir == 'Y') = dy(welldir == 'Y');
thicknessWell(welldir == 'X') = dx(welldir == 'X');

if ~isfield(W, 'rR')
    error('The representative radius of the well is not initialized');
end
rR = vertcat(W.rR);

velocW = Tw.*(pW - pBH)./(2 * pi * rR .* thicknessWell);

end


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
