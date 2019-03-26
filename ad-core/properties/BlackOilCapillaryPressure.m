classdef BlackOilCapillaryPressure < GridProperty
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
            prop = prop@GridProperty(model, varargin{:});
            prop = prop.dependsOn('s', 'state');
        end
        
        function pc = evaluateOnDomain(prop, model, state)
            [act, phInd] = model.getActivePhases();
            nph = sum(act);
            pc = cell(1, nph);
            
            f = model.fluid;
            if model.water && model.oil && isfield(f, 'pcOW')
                sW = model.getProp(state, 'sw');
                pcow = prop.evaluateFunctionOnGrid(f.pcOW, sW);
                if ~isempty(prop.endpointOptionSW)
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
                pc{phInd == 3} = prop.evaluateFunctionOnGrid(f.pcOG, sG);
            end
        end
        
        function anyPresent = pcPresent(prop, model)
            f = model.fluid;
            anyPresent = isfield(f, 'pcOW') || isfield(f, 'pcOG');
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
                pc = model.fluid.pcOW{i};
                v = pc([0; 1]);
                prop.pcOWMin(i) = min(v);
                prop.pcOWMax(i) = max(v);
            end
            prop.endpointOptionSW = option;
            pc_min = prop.evaluateFunctionOnGrid(model.fluid.pcOW, sw_prescribed);
            prop.endpointSW = sw_prescribed;
            prop.endpointPCOW = pc_min;
        end
    end
end