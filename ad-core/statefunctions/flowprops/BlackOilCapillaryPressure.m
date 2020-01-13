classdef BlackOilCapillaryPressure < StateFunction & SaturationProperty
    % Implementation of black-oil type capillary pressure
    properties

    end

    properties (Access = protected)

    end
    
    methods
        function prop = BlackOilCapillaryPressure(model, varargin)
            prop = prop@StateFunction(model, varargin{:});
            prop = prop.dependsOn('s', 'state');
        end
        
        function pc = evaluateOnDomain(prop, model, state)
            [act, phInd] = model.getActivePhases();
            nph = sum(act);
            pc = cell(1, nph);
            
            f = model.fluid;
            if model.water && model.oil && isfield(f, 'pcOW')
                sW = model.getProp(state, 'sw');
                if prop.scalingActive
%                     swcon = prop.getConnateWater(model, state);
                    pts = model.rock.krscale.drainage;
                    reg = prop.regions;
                    [get, ~, U, L] = SaturationProperty.getSatPointPicker(f, pts, reg, prop.cell_subset);
                    [swcon, SWCON] = get('w', L);
                    [swmax, SWMAX] = get('w', U);
                    sW = swcon + (sW - SWCON).*(swmax - swcon)./(SWMAX - SWCON);
                end
                pcow = prop.evaluateFunctionOnDomainWithArguments(f.pcOW, sW);
                if isfield(state, 'pcowScale')
                    pcow = pcow.*state.pcowScale;
                end
                % Note sign! Water is always first
                pc{phInd == 1} = -pcow;
            end
            
            if model.gas && model.oil && isfield(f, 'pcOG')
                sG = model.getProp(state, 'sg');
                pc{phInd == 3} = prop.evaluateFunctionOnDomainWithArguments(f.pcOG, sG);
            end
            if ~model.oil && isfield(f, 'pcWG')
                pc{phInd == 2} = prop.evaluateFunctionOnDomainWithArguments(f.pcWG, sG);
            end
        end
        
        function anyPresent = pcPresent(prop, model)
            f = model.fluid;
            anyPresent = isfield(f, 'pcOW') || isfield(f, 'pcOG') || isfield(f, 'pcWG');
        end
        
        function property = subset(property, subs)
            property = subset@StateFunction(property, subs);
            property.cell_subset = subs;
        end

    end
end

%{
Copyright 2009-2019 SINTEF Digital, Mathematics & Cybernetics.

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
