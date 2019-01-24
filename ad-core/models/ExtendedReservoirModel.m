classdef ExtendedReservoirModel
    properties
        
    end
    
    methods
        
        function c = getComponentDensity(model, state, name, ph)
            if nargin < 4
                ph = model.getPhaseNames();
            end
            nph = numel(ph);

            c = cell(nph, 1);
            switch name
                case 'water'
                    wix = strcmpi(ph, 'w');
                    c{wix} = 1;
                case 'oil'
                    oix = strcmpi(ph, 'o');
                    c{oix} = 1;
                case 'gas'
                    gix = strcmpi(ph, 'g');
                    c{gix} = 1;
                otherwise
                    error('Unknown component %s', name);
            end
        end
        
        function c = getPhaseComposition(model, state, name, varargin)
            c = model.getComponentDensity(state, name, varargin{:});
            rho = model.getProp(state, 'Density');
            c = cellfun(@rdivide, c, rho, 'UniformOutput', false);
        end
        
        function mass = getComponentMass(model, state, name, varargin)
            c = model.getComponentDensity(state, name, varargin{:});
            pv = model.getProp(state, 'PoreVolume');
            mass = cellfun(@(x) x.*pv, c, 'UniformOutput', false);
        end
        
        function mass = getComponentMobility(model, state, name, varargin)
            c = model.getComponentDensity(state, name, varargin{:});
            mob = model.getProp(state, 'Mobility');
            mass = cellfun(@(x) x.*mob, c, 'UniformOutput', false);
        end
    end
end