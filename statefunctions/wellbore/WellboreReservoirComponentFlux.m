classdef WellboreReservoirComponentFlux < CouplingTerm
%State function for mass flux between wellbore and reservoir
    
    properties
    end
    
    methods
        
        function prop = WellboreReservoirComponentFlux(model, reservoir, wellbore, varargin)
            
            % Optional input arguments
            opt = struct( ...
                'cellsReservoir', [], ...
                'cellsWellbore' , []  ...
            );
            [opt, extra] = merge_options(opt, varargin{:});
            % Parent class constructor
            prop@CouplingTerm(model, extra{:});
            
            % Check that models share the same component names
            assert(isa(model.submodels.(wellbore), 'WellboreModel'), ...
                'Wellbore model must be of class `WellboreModel`');
            assert( ...
                isa(model.submodels.(reservoir), 'GeothermalModel') && ...
                isa(model.submodels.(wellbore).parentModel, 'GeothermalModel'), ...
                ['State function currently only implemented for ', ...
                'Geothermal models']                               ...
            );
            
            % Set reservoir and wellbore cells if not given
            if isempty(opt.cellsReservoir)
                opt.cellsReservoir ...
                    = model.submodels.(wellbore).getGlobalWellCells();
            end
            if isempty(opt.cellsWellbore)
                opt.cellsWellbore ...
                    = model.submodels.(wellbore).G.cells.type == 0;
            end
            
            if isempty(prop.couplings)
                % Specify coupling
                coupling = { ...
                    struct( ...
                        'model'    , reservoir         , ... % Wellbore model name
                        'equations', {{'H2O'}}         , ... % Equations
                        'subset'   , opt.cellsReservoir, ... % Perforated reservoir cells
                        'sign'     , -1                 ... % Source
                    ), ...
                    struct( ...
                        'model'    , wellbore          , ... % Model name
                        'equations', {{'H2O'}}         , ... % Equations
                        'subset'   , opt.cellsWellbore , ... % Perforated well cells
                        'sign'     , 1                   ... % Sink
                    ), ...
                };
                % Set to state function
                prop.couplings = coupling;
                
            end
            prop.label = 'Q^{w,r}';
            prop.submodels.reservoir = reservoir;
            prop.submodels.well      = wellbore;
            
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

%{
Copyright 2009-2024 SINTEF Digital, Mathematics & Cybernetics.

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