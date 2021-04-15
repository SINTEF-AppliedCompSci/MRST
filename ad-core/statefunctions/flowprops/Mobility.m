classdef Mobility < StateFunction
    % Mobility for each phase. Normally rel. perm. divided by viscosity.
    properties (Access = protected)
        hasMultipliers = false;
    end
    
    methods
        function gp = Mobility(model, varargin)
            gp@StateFunction(model, varargin{:});
            gp = gp.dependsOn({'RelativePermeability'});
            gp = gp.dependsOn('Viscosity', 'PVTPropertyFunctions');
            gp.label = '\lambda_\alpha';
            gp.outputRange = [0, inf];
        end
        function mob = evaluateOnDomain(prop, model, state)
            kr = prop.getEvaluatedDependencies(state, 'RelativePermeability');
            mu = prop.getEvaluatedExternals(model, state, 'Viscosity');
            mob = cellfun(@(x, y) x./y, kr, mu, 'UniformOutput', false);
            if isfield(model.fluid, 'tranMultR')
                % Pressure dependent mobility multiplier 
                p = model.getProp(state, 'pressure');
                mult = model.fluid.tranMultR(p);
                mob = cellfun(@(x) x.*mult, mob, 'UniformOutput', false);
            end
            if prop.hasMultipliers
                mult = prop.getEvaluatedDependencies(state, 'MobilityMultipliers');
                for i = 1:numel(mult)
                    m = mult{i};
                    if ~isempty(m)
                        mu{i} = mu{i}.*m;
                    end
                end
            end
            % Check for negative values
            mv = cellfun(@(x) min(value(x)), mob);
            isAD = cellfun(@(x) isa(x, 'ADI'), mob);
            if any(isAD) && ~all(isAD)
                s = mob(isAD);
                s = s{1};
                bad = find(~isAD);
                for i = 1:numel(bad)
                    ix = bad(i);
                    mob{ix} = model.AutoDiffBackend.convertToAD(mob{ix}, s);
                end
            end
        end
        
        function prop = enableMultipliers(prop)
            prop.hasMultipliers = true;
            prop = prop.dependsOn('MobilityMultipliers');
        end
    end
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
