classdef WellboreReservoirComponentFlux < CouplingTerm
   
    properties
    end
    
    methods
        
        function prop = WellboreReservoirComponentFlux(model, well, reservoir, varargin)
            
            prop@CouplingTerm(model, varargin{:});
            prop.label = 'q';
            prop.submodels.well      = well;
            prop.submodels.reservoir = reservoir;
            
        end
        
        function Q = evaluateOnDomain(prop, model, state)
            
            reservoir = prop.submodels.reservoir;
            well      = prop.submodels.well;
            cells     = model.submodels.(well).getGlobalWellCells();
            isPerf    = model.submodels.(well).G.cells.type == 0;
            
            % Get reservoir pressure in perforated cells
            p = model.submodels.(reservoir).getProp(state.(reservoir), 'Pressure');
            p = p(cells);
            % Get well pressure
            pw = model.submodels.(well).getProp(state.(well), 'Pressure');
            pw = pw(isPerf);
            % Get well indices
            WI = model.submodels.(well).parentModel.operators.WI;
            % Compute well total flux
            q = -WI.*(p - pw);
            % Compute upwind flag
            flag = q >= 0;

            % Get reservoir mobility in perforated cells
            mob = model.submodels.(reservoir).getProp(state.(reservoir), 'ComponentMobility');
            mob(~cellfun(@isempty, mob)) = applyFunction(@(x) x(cells), mob(~cellfun(@isempty, mob)));
            % Get well pressure
            mobw = model.submodels.(well).parentModel.getProp(state.(well), 'ComponentMobility');
            mobw(~cellfun(@isempty, mobw)) = applyFunction(@(mob) mob(isPerf), mobw(~cellfun(@isempty, mobw)));
            
            % Compute component fluxes in to/out of wellbore out of/in to
            % reservoir (ComponentMobility is multiplied by density)
            ncomp = model.submodels.(reservoir).getNumberOfComponents();
            nph   = model.submodels.(reservoir).getNumberOfPhases();
            Q = cell(ncomp, 1);
            for c = 1:ncomp
                Q{c} = 0;
                for ph = 1:nph
                    if isempty(mob{c,ph}), continue; end
                    mobPh = flag.*mobw{c, ph}  + ~flag.*mob{c, ph};
                    Q{c} = Q{c} + mobPh.*q;
                end
            end
            
        end
            
    end
    
end