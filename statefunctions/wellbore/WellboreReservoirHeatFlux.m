classdef WellboreReservoirHeatFlux < CouplingTerm
   
    properties
    end
    
    methods
        
        function prop = WellboreReservoirHeatFlux(model, well, reservoir, varargin)
            
            prop@CouplingTerm(model, varargin{:});
            prop.label = 'q';
            prop.submodels.well      = well;
            prop.submodels.reservoir = reservoir;
            
        end
        
        function Qh = evaluateOnDomain(prop, model, state)
            
            Q = model.getProp(state, 'WellboreReservoirComponentFlux');
            
            reservoir = prop.submodels.reservoir;
            well      = prop.submodels.well;
            cells     = model.submodels.(well).getWellCells();
            isPerf    = model.submodels.(well).G.cells.isPerf;
            
            h = model.submodels.(reservoir).getProp(state.(reservoir), 'PhaseEnthalpy');
            h = h{1}(cells);
            
            hw = model.submodels.(well).getProp(state.(well), 'PhaseEnthalpy');
            hw = hw{1}; hw = hw(isPerf);
            flag = Q{1} >= 0;
            Qh = Q{1}.*(flag.*hw + ~flag.*h);
            
            % Get reservoir temperature in perforated cells
            T = model.submodels.(reservoir).getProp(state.(reservoir), 'Temperature');
            T = T(cells);
            % Get well temperature
            Tw = model.submodels.(well).getProp(state.(well), 'Temperature');
            Tw = Tw(isPerf);
            % Get well indices
            WIth = model.submodels.(well).parentModel.operators.WIth;
            % Compute well total flux
            Qh = Qh - WIth.*(T - Tw);
            
        end
            
    end
    
end