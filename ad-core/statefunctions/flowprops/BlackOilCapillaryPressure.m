classdef BlackOilCapillaryPressure < StateFunction
    properties
    end
    
    properties (Access = protected)
        endpointSW
        endpointPCOW
        pcOWMin
        pcOWMax
        endpointOptionSW = [];
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
                pcow = prop.evaluateFunctionOnDomainWithArguments(f.pcOW, sW);
                if isfield(state, 'pcowScale')
                    pcow = pcow.*state.pcowScale;
                end
                if ~isempty(prop.endpointOptionSW) && false % Disabled due to bad data
                    pcmin = prop.pcOWMin;
                    pcmax = prop.pcOWMax;
                    pcw = prop.endpointPCOW;
                    reg = prop.regions;
                    if ~isempty(reg)
                        pcmin = pcmin(reg);
                        pcmax = pcmax(reg);
                        pcw = pcw(reg);
                    end
                    pc_scale = (pcow - pcmin)./(pcmax - pcmin);
                    switch prop.endpointOptionSW
                        case 1
                            % Initial water is interpreted as maximum pc
                            pcow = (pcw - pcmin).*pc_scale + pcmin;
                        case 2
                            pcow = (pcmax - pcw).*pc_scale + pcw;
                    end
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
        
        function prop = setWaterEndpointScaling(prop, model, sw_prescribed, option)
            % Special case where water capillary pressure is adjusted to
            % match the initial water saturation, instead of the other way
            % around
            if ~isfield(model.fluid, 'pcOW')
                return
            end
            nreg = numel(model.fluid.pcOW);
            prop.endpointSW = zeros(nreg, 1);
            prop.endpointPCOW = zeros(nreg, 1);
            prop.pcOWMin = zeros(nreg, 1);
            prop.pcOWMax = zeros(nreg, 1);
            for i = 1:nreg
                if isa(model.fluid.pcOW, 'function_handle')
                    pc = model.fluid.pcOW;
                else
                    pc = model.fluid.pcOW{i};
                end
                v = pc([0; 1]);
                prop.pcOWMin(i) = min(v);
                prop.pcOWMax(i) = max(v);
            end
            prop.endpointOptionSW = option;
            pc_min = prop.evaluateFunctionOnDomainWithArguments(model.fluid.pcOW, sw_prescribed);
            prop.endpointSW = sw_prescribed;
            prop.endpointPCOW = pc_min;
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
