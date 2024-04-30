classdef WellboreReservoirHeatFlux < CouplingTerm
%State function for heat flux between wellbore and reservoir

    properties
    end
    
    methods
        
        function prop = WellboreReservoirHeatFlux(model, reservoir, wellbore, varargin)
            
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
                        'equations', {{'energy'}}      , ... % Equations
                        'subset'   , opt.cellsReservoir, ... % Perforated reservoir cells
                        'sign'     , -1                 ... % Source
                    ), ...
                    struct( ...
                        'model'    , wellbore          , ... % Model name
                        'equations', {{'energy'}}      , ... % Equations
                        'subset'   , opt.cellsWellbore , ... % Perforated well cells
                        'sign'     , 1                   ... % Sink
                    ), ...
                };
                % Set to state function
                prop.couplings = coupling;
                
            end
            prop.label = 'Q_h^{w,r}';
            prop.submodels.reservoir = reservoir;
            prop.submodels.well      = wellbore;
            
        end
        
        function Qh = evaluateOnDomain(prop, model, state)
            
            % Get mass flux between well and reservoir
            Q = model.getProp(state, 'WellboreReservoirComponentFlux');
            
            reservoir = prop.submodels.reservoir;
            well      = prop.submodels.well;
            cells     = model.submodels.(well).getGlobalWellCells();
            isPerf    = model.submodels.(well).G.cells.type == 0;
            
            % Compute enthalpy in perforated cells of reservoir
            h = model.submodels.(reservoir).getProp(state.(reservoir), 'PhaseEnthalpy');
            h = h{1}(cells);
            % Compute enthalpy in perforated segments of wellbore
            hw = model.submodels.(well).getProp(state.(well), 'PhaseEnthalpy');
            hw = hw{1}; hw = hw(isPerf);
            % Upstream weight enthalpy and compute advective heat flux
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
            % Compute well total flux as sum af advective and conductive
            % heat flux
            Qh = Qh - WIth.*(T - Tw);
            
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