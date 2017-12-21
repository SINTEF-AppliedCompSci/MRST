classdef OilWaterSurfactantModel < TwoPhaseOilWaterModel
    %
    %
    % SYNOPSIS:
    %   model = FullyImplicitOilWaterSurfactantModel(G, rock, fluid, varargin)
    %
    % DESCRIPTION: 
    %   Fully implicit model for an oil water system with surfactant. All
    %   the equations are solved implicitly. A description of the surfactant model
    %   that is implemented can be found in the directory ad-eor/docs .
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
    % SEE ALSO: equationsOilWaterSurfactant, ImplicitExplicitOilWaterSurfactantModel
    %

    properties
        surfactant
    end

    methods

        function model = OilWaterSurfactantModel(G, rock, fluid, varargin)

            model = model@TwoPhaseOilWaterModel(G, rock, fluid, varargin{:});
            model = model.setupOperators(G, rock, varargin{:});
            model.surfactant = true;
            model = merge_options(model, varargin{:});

        end

        function model = setupOperators(model, G, rock, varargin)
            model.operators.veloc = computeVelocTPFA(G, model.operators.internalConn);
            model.operators.sqVeloc = computeSqVelocTPFA(G, model.operators.internalConn);
        end

        function [problem, state] = getEquations(model, state0, state, dt, drivingForces, varargin)
            [problem, state] = equationsOilWaterSurfactant(state0, state, model, dt, drivingForces, ...
                                                           varargin{:});
        end

        function [state, report] = updateState(model, state, problem, dx, drivingForces)
            [state, report] = updateState@TwoPhaseOilWaterModel(model, state, problem,  dx, ...
                                                              drivingForces);
            % cap the concentration (only if implicit solver for concentration)
            if model.surfactant
                c = model.getProp(state, 'surfactant');
                state = model.setProp(state, 'surfactant', max(c, 0) );
            end
        end

        function [state, report] = updateAfterConvergence(model, state0, state, dt, drivingForces)
            [state, report] = updateAfterConvergence@TwoPhaseOilWaterModel(model, state0, state, dt, ...
                                                              drivingForces);
              if model.surfactant
                  c     = model.getProp(state, 'surfactant');
                  cmax  = model.getProp(state, 'surfactantmax');
                  state = model.setProp(state, 'surfactantmax', max(cmax, c));
              end
        end

        function varargout = evaluateRelPerm(model, sat, varargin)
            error('function evaluateRelPerm is not implemented for surfactant model')
        end

        function state = validateState(model, state)
            state = validateState@TwoPhaseOilWaterModel(model, state);
            nc = model.G.cells.num;
            model.checkProperty(state, 'Surfactant', [nc, 1], [1, 2]);
            model.checkProperty(state, 'SurfactantMax', [nc, 1], [1, 2]);
        end

        function [fn, index] = getVariableField(model, name)
            switch(lower(name))
                case {'surfactant'}
                    index = 1;
                    fn = 'c';
                case {'surfactantmax'}
                    index = 1;
                    fn = 'cmax';
                case 'qwsft'
                    index = 1;
                    fn = 'qWSft';
                otherwise
                    [fn, index] = getVariableField@TwoPhaseOilWaterModel(...
                        model, name);
            end
        end

        function names = getComponentNames(model)
            names = getComponentNames@TwoPhaseOilWaterModel(model);
            if model.surfactant
                names{end+1} = 'surfactant';
            end
        end


        function state = storeSurfData(model, state, s, c, Nc, sigma)
            state.SWAT    = double(s);
            state.SURFACT = double(c);
            state.SURFCNM = log(double(Nc))/log(10);
            state.SURFST  = double(sigma);
            % state.SURFADS = double(ads);
        end

        function [names, types] = getExtraWellEquationNames(model)
            [names, types] = getExtraWellEquationNames@TwoPhaseOilWaterModel(model);
            if model.surfactant
                names{end+1} = 'surfactantWells';
                types{end+1} = 'perf';
            end
        end

        function names = getExtraWellPrimaryVariableNames(model)
            names = getExtraWellPrimaryVariableNames@TwoPhaseOilWaterModel(model);
            if model.surfactant
                names{end+1} = 'qWSft';
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
              case {'surfactant'}
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
        
        function [compEqs, compSrc, compNames, wellSol] = getExtraWellContributions(model, well, wellSol0, wellSol, q_s, bh, packed, qMass, qVol, dt, iteration)
            [compEqs, compSrc, compNames, wellSol] = getExtraWellContributions@TwoPhaseOilWaterModel(model, well, wellSol0, wellSol, q_s, bh, packed, qMass, qVol, dt, iteration);
            if model.surfactant
                % Implementation of surfactant source terms.
                %
                assert(model.water, 'Surfactant injection requires a water phase.');
                f = model.fluid;
                if well.isInjector
                    concWell = model.getProp(well.W, 'surfactant');
                else
                    pix = strcmpi(model.getComponentNames(), 'surfactant');
                    concWell = packed.components{pix};
                end
                qwsft = packed.extravars{strcmpi(packed.extravars_names, 'qwsft')};
                % Water is always first
                wix = 1;
                cqWs = qMass{wix}./f.rhoWS; % get volume rate, at
                                            % surface condition.
                cqS = concWell.*cqWs;

                compEqs{end+1} = qwsft - sum(cqWs);
                compSrc{end+1} = cqS;
                compNames{end+1} = 'surfactantWells';
            end
        end
    end
end


%{
Copyright 2009-2017 SINTEF ICT, Applied Mathematics.

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

