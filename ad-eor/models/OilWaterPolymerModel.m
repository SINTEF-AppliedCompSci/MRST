classdef OilWaterPolymerModel < TwoPhaseOilWaterModel
    % Oil/water/polymer system
    %
    %
    % SYNOPSIS:
    %   model = OilWaterPolymerModel(G, rock, fluid, varargin)
    %
    % DESCRIPTION:
    %   Two phase model with polymer. A description of the polymer model
    %   that is implemented here can be found in the directory ad-eor/docs .
    %
    % PARAMETERS:
    %   G        - Grid
    %   rock     - Rock structure
    %   fluid    - Fluid structure
    %   varargin - optional parameter
    %
    % RETURNS:
    %   class instance
    %
    % EXAMPLE:
    %
    % SEE ALSO: ThreePhaseBlackOilPolymerModel
    %

    properties
        % Polymer present
        polymer

    end

    methods
        function model = OilWaterPolymerModel(G, rock, fluid, varargin)

            model = model@TwoPhaseOilWaterModel(G, rock, fluid);
            % This is the model parameters for oil/water/polymer
            model.polymer = true;
            model = merge_options(model, varargin{:});

        end

        function [problem, state] = getEquations(model, state0, state, ...
                dt, drivingForces, varargin)
            [problem, state] = equationsOilWaterPolymer(state0, state, ...
                model, dt, drivingForces, varargin{:});
        end

        function state = validateState(model, state)
            state = validateState@TwoPhaseOilWaterModel(model, state);
            % Polymer must be present
            model.checkProperty(state, 'Polymer', model.G.cells.num, 1);
            fn = model.getVariableField('polymermax');
            if ~isfield(state, fn)
                state.(fn) = model.getProp(state, 'Polymer');
            end
        end

        function [state, report] = updateState(model, state, problem, ...
                dx, drivingForces)

            if model.polymer
                % Store the polymer from previous iteration temporarily to
                % use in convergence criteria
                cp_prev = model.getProp(state, 'polymer');
            end

            [state, report] = updateState@TwoPhaseOilWaterModel(model, ...
               state, problem,  dx, drivingForces);

            if model.polymer
                % Limit polymer concentration to [0, fluid.cpmax]
                cp = model.getProp(state, 'polymer');
                cp = min(cp, model.fluid.cpmax);
                state = model.setProp(state, 'polymer', max(cp, 0) );
                state.cp_prev = cp_prev;

                % Shear Thinning Report
                % We (may) have stored the shear thinning report
                % temporarily in the state strucpture. We move this over to
                % the report structure instead. The reason for this is that
                % there is no report returned from the equations.
                if isfield(state, 'ShearThinningReport')
                    report.ShearThinning = state.ShearThinningReport;
                    state = rmfield(state, 'ShearThinningReport');
                end
            end
        end

        function [state, report] = updateAfterConvergence(model, state0, state, dt, drivingForces)
            [state, report] = updateAfterConvergence@TwoPhaseOilWaterModel(model, state0, state, dt, drivingForces);
            if model.polymer
                cp    = model.getProp(state, 'polymer');
                cpmax = model.getProp(state, 'polymermax');
                state = model.setProp(state, 'polymermax', max(cpmax, cp));

                if isfield(state, 'cp_prev')
                    % Remove the temporary field used for convergence
                    state = rmfield(state, 'cp_prev');
                end
            end
        end


        function [fn, index] = getVariableField(model, name, varargin)
            % Get the index/name mapping for the model (such as where
            % pressure or water saturation is located in state)
            index = 1;
            switch(lower(name))
                case {'polymer', 'polymermax'}
                    c = model.getComponentNames();
                    index = find(strcmpi(c, 'polymer'));
                    if strcmpi(name, 'polymer')
                        fn = 'cp';
                    else
                        fn = 'cpmax';
                    end
                otherwise
                    [fn, index] = getVariableField@TwoPhaseOilWaterModel(...
                                    model, name, varargin{:});
            end
        end

        function names = getComponentNames(model)
            names = getComponentNames@TwoPhaseOilWaterModel(model);
            if model.polymer
                names{end+1} = 'polymer';
            end
        end


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

        function [compEqs, compSrc, eqNames, wellSol] = getExtraWellContributions(model, well, wellSol0, wellSol, q_s, bh, packed, qMass, qVol, dt, iteration)
            [compEqs, compSrc, eqNames, wellSol] = getExtraWellContributions@TwoPhaseOilWaterModel(model, well, wellSol0, wellSol, q_s, bh, packed, qMass, qVol, dt, iteration);

            if model.polymer
                assert(model.water, 'Polymer injection requires a water phase.');
                f = model.fluid;

                % Polymer concentration given by the well control
                concpolyCtrl = model.getProp(well.W, 'polymer');

                % Polymer concentration from the reservoir at the connections
                pix = strcmpi(model.getComponentNames(), 'polymer');
                concpolyRes = packed.components{pix};

                % Water is always first
                wix = 1;
                cqWs = qMass{wix}./f.rhoWS; % connection volume flux at surface condition
                q_sw = q_s{wix}; % water flow from bottom hole. 

                isInj = (cqWs > 0);
                % Computation of total polymer influx in well qP
                % (We recall that cqWs<0 when producing)
                qP = sum(-concpolyRes(~isInj).*cqWs(~isInj));
                if (q_sw > 0) 
                    % Well is injecting water
                    qP = qP + concpolyCtrl*q_sw;
                end
                
                % Computation of total water outflux from well
                qWout = sum(cqWs(isInj));
                if (q_sw < 0) 
                    % Well is injecting water
                    % (We recall that q_sw < 0 when producing)
                    qWout = qWout - q_sw;
                end

                % Let us now compute the polymer flow rate at each connection
                %
                % First, we handle the producting connections.
                % The term (a + (1 - a).*cbarw) account for the todd-longstaff mixing factor,
                % which model the fact that for not-fully mixed polymer solution
                % the polymer does not travel at the same velocity as water. See
                % the governing equation for polymer
                % (e.g. equationsOilWaterPolymer.m)
                cbarw = concpolyRes/f.cpmax;
                a     = f.muWMult(f.cpmax).^(1-f.mixPar);
                cqP   = concpolyRes.*cqWs./(a + (1-a).*cbarw);

                % For the injecting connections, the polymer is distributed equally in each
                % connection.

                if any(isInj)
                    cqP(isInj) = (qP./qWout).*cqWs(isInj);
                end                
                
                compSrc{end+1} = cqP;
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
