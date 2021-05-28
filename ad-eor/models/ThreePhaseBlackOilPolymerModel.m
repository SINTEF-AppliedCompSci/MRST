classdef ThreePhaseBlackOilPolymerModel < ThreePhaseBlackOilModel
    %Three-phase black-oil model with support for polymer injection
    %
    % SYNOPSIS:
    %   model = ThreePhaseBlackOilPolymerModel(G, rock, fluid, varargin)
    %
    % DESCRIPTION: 
    %   Fully implicit three phase blackoil model with polymer.
    %
    % PARAMETERS:
    %   G        - Grid
    %   rock     - Rock structure
    %   fluid    - Fluid structure
    %   varargin - optional parameters
    %
    % RETURNS:
    %   class instance
    %
    % EXAMPLE:
    %
    % SEE ALSO:  equationsThreePhaseBlackOilPolymer, OilWaterPolymerModel
    %

    properties
        % Polymer present
        polymer
        % Using PLYSHEAR shear model based on water velocity
        usingShear
	    % Using PLYSHLOG shear model based on water velocity
        usingShearLog
	    % Using PLYSHLOG shear model base on water shear rate
        usingShearLogshrate
    end

    methods
        function model = ThreePhaseBlackOilPolymerModel(G, rock, fluid, varargin)
            model = model@ThreePhaseBlackOilModel(G, rock, fluid, varargin{:});

            % This is the model parameters for oil/water/gas/polymer system
            model.polymer = true;

            if (isfield(fluid,'shrate') && ~isfield(fluid,'plyshlog'))
                error('SHRATE is specified while PLYSHLOG is not specified')
            end
            if (isfield(fluid, 'plyshearMult') && isfield(fluid, 'plyshlog'))
                error('PLYSHLOG and PLYSHEAR are existing together');
            end

            model.usingShear = isfield(fluid, 'plyshearMult');
            model.usingShearLog = (isfield(fluid, 'plyshlog') && ~isfield(fluid, 'shrate'));
            model.usingShearLogshrate = (isfield(fluid, 'plyshlog') && isfield(fluid, 'shrate'));
            % model.wellVarNames = {'qWs', 'qOs', 'qGs', 'qWPoly', 'bhp'};
            model = merge_options(model, varargin{:});
        end

        % --------------------------------------------------------------------%
        function [problem, state] = getEquations(model, state0, state, dt, drivingForces, varargin)
            [problem, state] = equationsThreePhaseBlackOilPolymer(state0, state, ...
                model, dt, drivingForces, varargin{:});
        end

        % --------------------------------------------------------------------%
        function state = validateState(model, state)
            state = validateState@ThreePhaseBlackOilModel(model, state);
            % Polymer must be present
            model.checkProperty(state, 'Polymer', [model.G.cells.num, 1], [1, 2]);
        end

        % --------------------------------------------------------------------%
        function [state, report] = updateState(model, state, problem, dx, drivingForces)
            [state, report] = updateState@ThreePhaseBlackOilModel(model, ...
               state, problem,  dx, drivingForces);

            if model.polymer
                cp = model.getProp(state, 'polymer');
                cp = min(cp, model.fluid.cpmax);
                state = model.setProp(state, 'polymer', max(cp, 0));
            end
        end

        % --------------------------------------------------------------------%
        function [state, report] = updateAfterConvergence(model, state0, state, dt, drivingForces)
            [state, report] = updateAfterConvergence@ThreePhaseBlackOilModel(model, state0, state, dt, drivingForces);
            if model.polymer
                cp     = model.getProp(state, 'polymer');
                cpmax  = model.getProp(state, 'polymermax');
                state = model.setProp(state, 'polymermax', max(cpmax, cp));
            end
        end

        % --------------------------------------------------------------------%
        function [fn, index] = getVariableField(model, name, varargin)
            % Get the index/name mapping for the model (such as where
            % pressure or water saturation is located in state)
            switch(lower(name))
                case {'polymer', 'polymermax'}
                    c = model.getComponentNames();
                    index = find(strcmpi(c, 'polymer'));
                    if strcmpi(name, 'polymer')
                        fn = 'cp';
                    else
                        fn = 'cpmax';
                    end
                case 'qwpoly'
                    fn = 'qWPoly';
                    index = ':';
                otherwise
                    [fn, index] = getVariableField@ThreePhaseBlackOilModel(...
                                    model, name, varargin{:});
            end
        end

        % --------------------------------------------------------------------%
        function names = getComponentNames(model)
            names = getComponentNames@ThreePhaseBlackOilModel(model);
            if model.polymer
                names{end+1} = 'polymer';
            end
        end
        % --------------------------------------------------------------------%
        function scaling = getScalingFactorsCPR(model, problem, names, solver)
            nNames = numel(names);

            scaling = cell(nNames, 1);
            handled = false(nNames, 1);

            for iter = 1:nNames
                name = lower(names{iter});
                switch name
                    case 'polymer'
                        s = 0;
                    otherwise
                        continue
                end
                sub = strcmpi(problem.equationNames, name);

                scaling{iter} = s;
                handled(sub) = true;
            end
            if ~all(handled)
                % Get rest of scaling factors
                other = getScalingFactorsCPR@ThreePhaseBlackOilModel(model, problem, names(~handled), solver);
                [scaling{~handled}] = other{:};
            end
        end

        function [eq, src] = addComponentContributions(model, cname, eq, component, src, force)
        % For a given component conservation equation, compute and add in
        % source terms for a specific source/bc where the fluxes have
        % already been computed.
        %
        % PARAMETERS:
        %
        %   model  - (Base class, automatic)
        %
        %   cname  - Name of the component. Must be a property known to the
        %            model itself through `getProp` and `getVariableField`.
        %
        %   eq     - Equation where the source terms are to be added. Should
        %            be one value per cell in the simulation grid (model.G)
        %            so that the src.sourceCells is meaningful.
        %
        %   component - Cell-wise values of the component in question. Used
        %               for outflow source terms only.
        %
        %   src    - Source struct containing fields for fluxes etc. Should
        %            be constructed from force and the current reservoir
        %            state by `computeSourcesAndBoundaryConditionsAD`.
        %
        %   force  - Force struct used to produce src. Should contain the
        %            field defining the component in question, so that the
        %            inflow of the component through the boundary condition
        %            or source terms can accurately by estimated.
            if isempty(force)
                return
            end
            c = model.getProp(force, cname);
            cells = src.sourceCells;
            switch lower(cname)
              case {'polymer'}
                % Water based EOR, multiply by water flux divided by
                % density and add into corresponding equation
                qW = src.phaseMass{1}./model.fluid.rhoWS;
                isInj = qW > 0;
                qC = (isInj.*c + ~isInj.*component(cells)).*qW;
              otherwise
                error(['Unknown component ''', cname, '''. BC not implemented.']);
            end
            eq(cells) = eq(cells) - qC;
            src.components{end+1} = qC;
        end
        
        function [names, types] = getExtraWellEquationNames(model)
            [names, types] = getExtraWellEquationNames@ThreePhaseBlackOilModel(model);
            if model.polymer
                names{end+1} = 'polymerWells';
                types{end+1} = 'perf';
            end
        end

        function names = getExtraWellPrimaryVariableNames(model)
            names = getExtraWellPrimaryVariableNames@ThreePhaseBlackOilModel(model);
            if model.polymer
                names{end+1} = 'qWPoly';
            end
        end

        function [compEqs, compSrc, eqNames, wellSol] = getExtraWellContributions(model, well, wellSol0, wellSol, q_s, bh, packed, qMass, qVol, dt, iteration)
            [compEqs, compSrc, eqNames, wellSol] = getExtraWellContributions@ThreePhaseBlackOilModel(model, well, wellSol0, wellSol, q_s, bh, packed, qMass, qVol, dt, iteration);
            if model.polymer
                assert(model.water, 'Polymer injection requires a water phase.');
                f = model.fluid;

                % Water is always first
                wix = 1;
                cqWs = qMass{wix}./f.rhoWS; % connection volume flux at surface condition

                if well.isInjector()
                    concWellp = model.getProp(well.W, 'polymer');
                    cqP = concWellp.*cqWs;
                else
                    pix = strcmpi(model.getComponentNames(), 'polymer');
                    concWellp = packed.components{pix};

                    a = f.muWMult(f.cpmax).^(1-f.mixPar);
                    cpbarw = concWellp/f.cpmax;

                    % the term (a + (1 - a).*cpbarw) account for the
                    % todd-longstaff mixing factor, which model the fact that for
                    % not-fully mixed polymer solution the polymer does not
                    % travel at the same velocity as water. See the governing
                    % equation for polymer (e.g. equationsOilWaterPolymer.m)
                    cqP = concWellp.*cqWs./(a + (1-a).*cpbarw);
                end

                qwpoly = packed.extravars{strcmpi(packed.extravars_names, 'qwpoly')};

                compEqs{end+1} = qwpoly - sum(concWellp.*cqWs);
                compSrc{end+1} = cqP;
                eqNames{end+1} = 'polymerWells';
            end
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