classdef WellConductiveHeatFlux < StateFunction
%State function for conductive heat flux between wellbore and reservoir

    properties
    end
    
    methods
        %-----------------------------------------------------------------%
        function gp = WellConductiveHeatFlux(model, varargin)
            gp@StateFunction(model, varargin{:});
            gp = gp.dependsOn({'FacilityWellMapping'});
            gp = gp.dependsOn({'temperature'}, 'state');
            gp = gp.dependsOn({'RockHeatTransmissibility', 'FluidHeatTransmissibility'}, 'FlowDiscretization');
            gp.label = 'q_{h,c}';
        end
        
        %-----------------------------------------------------------------%
        function qh = evaluateOnDomain(prop, model, state)
            
            map = prop.getEvaluatedDependencies(state, 'FacilityWellMapping');
            
            [T, lambdaR, lambdaF] = model.ReservoirModel.getProps(state, ...
                                            'Temperature'             , ...
                                            'RockThermalConductivity' , ...
                                            'FluidThermalConductivity');
            T   = T(map.cells);
            isT = strcmpi(state.FacilityState.names, 'well_temperature');
            if any(isT)
                Tw = state.FacilityState.primaryVariables{isT};
            else
                Tw = vertcat(map.W.T);
            end
            Tw = Tw(map.perf2well);
            if model.ReservoirModel.dynamicHeatTransFluid || ...
                    model.ReservoirModel.dynamicHeatTransRock
                % Assume thermal conductivity tensors are diagonal for now
                G      = model.ReservoirModel.G;
                lambda = lambdaR + lambdaF;
                rock   = struct('perm', ones(numel(value(lambda)), 1));
                wi     = cell(numel(map.W),1);
                for i = 1:numel(map.W)
                    radius = map.W(i).r;
                    cells  = map.W(i).cells;
                    if isfield(G.faces, 'nodePos')
                        % Compute well index
                        % Ensure equivalent greater than well radius
                        [dx, dy, dz] = cellDims(G, cells);
                        re = 2*0.14*sqrt(dx.^2 + dy.^2)/2;
                        C  = max(1.1*map.W(i).r(1)./re,1);
                        dx = [dx.*C, dy.*C, dz];
                    else
                        dx = [];
                    end
                    wi{i} = computeWellIndex(G, rock, radius, cells, 'cellDims', dx);
                    wi{i} = wi{i}.*lambda(cells);
                end
                wi = vertcat(wi{:});
            else
                wi = vertcat(map.W.WIth);
            end
            qh = -wi.*(T - Tw);
        end
    end
    
end

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