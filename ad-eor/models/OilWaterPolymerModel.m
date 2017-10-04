classdef OilWaterPolymerModel < TwoPhaseOilWaterModel
% Oil/water/polymer system
%
%
% SYNOPSIS:
%   model = OilWaterPolymerModel(G, rock, fluid, varargin)
%
% DESCRIPTION: Two phase model with polymer. A description of the polymer model
% that is implemented here can be found in the directory ad-eor/docs .
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

    properties
        % Polymer present
        polymer
        
        % Polymer differene tolerance
        useIncPolymerConvergence
        toleranceIncPolymer
        
        % Add extra output to wellsol/states for polymer quantities
        extraPolymerOutput
        
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
                c_prev = model.getProp(state, 'polymer');
            end
            
            [state, report] = updateState@TwoPhaseOilWaterModel(model, ...
               state, problem,  dx, drivingForces);

            if model.polymer
                % Limit polymer concentration to [0, fluid.cmax]
                c = model.getProp(state, 'polymer');
                c = min(c, model.fluid.cmax);
                state = model.setProp(state, 'polymer', max(c, 0) );
                state.c_prev = c_prev;
                
                % Shear Thinning Report               
                % We (may) have stored the shear thinning report
                % temporarily in the state structure. We move this over to
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
                c     = model.getProp(state, 'polymer');
                cmax  = model.getProp(state, 'polymermax');
                state = model.setProp(state, 'polymermax', max(cmax, c));
                
                if isfield(state, 'c_prev')
                    % Remove the temporary field used for convergence
                    state = rmfield(state, 'c_prev');
                end
            end
        end


        function [fn, index] = getVariableField(model, name)
            % Get the index/name mapping for the model (such as where
            % pressure or water saturation is located in state)
            index = 1;
            switch(lower(name))
                case {'polymer', 'polymermax'}
                    c = model.getComponentNames();
                    index = find(strcmpi(c, 'polymer'));
                    if strcmpi(name, 'polymer')
                        fn = 'c';
                    else
                        fn = 'cmax';
                    end
                case 'qwpoly'
                    fn = 'qWPoly';
                otherwise
                    [fn, index] = getVariableField@TwoPhaseOilWaterModel(...
                                    model, name);
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

        function [names, types] = getExtraWellEquationNames(model)
            [names, types] = getExtraWellEquationNames@TwoPhaseOilWaterModel(model);
            if model.polymer
                names{end+1} = 'polymerWells';
                types{end+1} = 'perf';
            end
        end

        function names = getExtraWellPrimaryVariableNames(model)
            names = getExtraWellPrimaryVariableNames@TwoPhaseOilWaterModel(model);
            if model.polymer
                names{end+1} = 'qWPoly';
            end
        end
        
        function [compEqs, compSrc, eqNames, wellSol] = getExtraWellContributions(model, well, wellSol0, wellSol, q_s, bh, packed, qMass, qVol, dt, iteration)
            [compEqs, compSrc, eqNames, wellSol] = getExtraWellContributions@TwoPhaseOilWaterModel(model, well, wellSol0, wellSol, q_s, bh, packed, qMass, qVol, dt, iteration);
            if model.polymer
                assert(model.water, 'Polymer injection requires a water phase.');
                f = model.fluid;
                if well.isInjector
                    concWell = model.getProp(well.W, 'polymer');
                else
                    pix = strcmpi(model.getComponentNames(), 'polymer');
                    concWell = packed.components{pix};
                end
                qwpoly = packed.extravars{strcmpi(packed.extravars_names, 'qwpoly')};
                a = f.muWMult(f.cmax).^(1-f.mixPar);
                cbarw     = concWell/f.cmax;
                % Water is always first
                wix = 1;
                cqWs = qMass{wix}./f.rhoWS; % connection volume flux at surface condition

                % the term (a + (1 - a).*cbarw) account for the
                % todd-longstaff mixing factor, which model the fact that for
                % not-fully mixed polymer solution the polymer does not
                % travel at the same velocity as water. See the governing
                % equation for polymer (e.g. equationsOilWaterPolymer.m)
                cqP = concWell.*cqWs./(a + (1-a).*cbarw);

                compEqs{end+1} = qwpoly - sum(concWell.*cqWs);
                compSrc{end+1} = cqP;
                eqNames{end+1} = 'polymerWells';
            end
        end
    end
end

