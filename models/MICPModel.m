classdef MICPModel < TwoPhaseOilWaterModel
% Script to implement the MICP model.
% 
% This script is modified from a file in The MATLAB Reservoir Simulation
% Toolbox (MRST), see
%   mrst/modules/ad-eor/models/OilWaterPolymerModel.m
%
% We refer to that script for a complete commented version of the file. 

%{ 
Partial copyright 2009-2021, SINTEF Digital, Mathematics & Cybernetics.
Partial copyright 2021, NORCE Norwegian Research Centre AS, Computational 
Geosciences and Modeling.

This file is part of the ad-micp module.

ad-micp is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.  

ad-micp is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this file.  If not, see <http://www.gnu.org/licenses/>.
%}
    properties
        % Microorganism present
        microorganism
        % Oxygen present
        oxygen
        % Urea present
        urea
        % Biofilm present
        biofilm
        % Calcite present
        calcite
    end

    methods
        function model = MICPModel(G, rock, fluid, varargin)
            model = model@TwoPhaseOilWaterModel(G, rock, fluid);
            % These are the model variables for MICP
            model.microorganism = true;
            model.oxygen = true;
            model.urea = true;
            model.biofilm = true;
            model.calcite = true;
            model = merge_options(model, varargin{:});
        end

        function [problem, state] = getEquations(model, state0, state, ...
                                               dt, drivingForces, varargin)
            [problem, state] = equationsMICP(state0, state, ...
                                    model, dt, drivingForces, varargin{:});
        end

        function state = validateState(model, state)
            state = validateState@TwoPhaseOilWaterModel(model, state);
            % Components must be present
            model.checkProperty(state, 'Microorganism', ...
                                                     model.G.cells.num, 1);
            model.checkProperty(state, 'Oxygen', model.G.cells.num, 1);
            model.checkProperty(state, 'Urea', model.G.cells.num, 1);
            model.checkProperty(state, 'Biofilm', model.G.cells.num, 1);
            model.checkProperty(state, 'Calcite', model.G.cells.num, 1);
        end

        function [state, report] = updateState(model, state, problem, ...
                dx, drivingForces)
            % Store the variables from previous iteration temporarily to
            % use in convergence criteria
            m_prev = model.getProp(state, 'microorganism');
            o_prev = model.getProp(state, 'oxygen');
            u_prev = model.getProp(state, 'urea');
            b_prev = model.getProp(state, 'biofilm');
            c_prev = model.getProp(state, 'calcite');

            [state, report] = updateState@TwoPhaseOilWaterModel(model, ...
                                       state, problem,  dx, drivingForces);
            % Limit the variables to [0, cmax]. For the microorganisms we 
            % set this value equal to the maximum biomass density value
            % we found in literature (105 kg/m^3). For the oxygen and urea
            % we set this value to the maximum injected concentration. For
            % the biofilm and calcite, we set this value equal to the
            % porosity minus the tolerance for the clogging criteria. 
            m = model.getProp(state, 'microorganism');
            m = min(m, 105);
            state = model.setProp(state, 'microorganism', max(m, 0) );
            state.m_prev = m_prev;
            o = model.getProp(state, 'oxygen');
            o = min(o, model.fluid.Comax);
            state = model.setProp(state, 'oxygen', max(o, 0) );
            state.o_prev = o_prev;
            u = model.getProp(state, 'urea');
            u = min(u, model.fluid.Cumax);
            state = model.setProp(state, 'urea', max(u, 0) );
            state.u_prev = u_prev;
            b = model.getProp(state, 'biofilm');
            b = min(b, model.rock.poro - model.fluid.ptol);
            state = model.setProp(state, 'biofilm', max(b, 0) );
            state.b_prev = b_prev;
            c = model.getProp(state, 'calcite');
            c = min(c, model.rock.poro - model.fluid.ptol);
            state = model.setProp(state, 'calcite', max(c, 0) );
            state.c_prev = c_prev;
        end

        function [state, report] = updateAfterConvergence(model, ...
                                          state0, state, dt, drivingForces)
            [state, report] = ...
            updateAfterConvergence@TwoPhaseOilWaterModel(model, state0, ...
                                                 state, dt, drivingForces);
                    % Remove the temporary field used for convergence
                    state = rmfield(state, 'm_prev');
                    state = rmfield(state, 'o_prev');
                    state = rmfield(state, 'u_prev');
                    state = rmfield(state, 'b_prev');
                    state = rmfield(state, 'c_prev');
        end

        function [fn, index] = getVariableField(model, name, varargin)
            % Get the index/name mapping for the model
            switch(lower(name))
                case {'microorganism'}
                    m = model.getComponentNames();
                    index = find(strcmpi(m, 'microorganism'));
                    fn = 'm';
                case {'oxygen'}
                    o = model.getComponentNames();
                    index = find(strcmpi(o, 'oxygen')) - 1;
                    fn = 'o';
                case {'urea'}
                    u = model.getComponentNames();
                    index = find(strcmpi(u, 'urea')) - 2;
                    fn = 'u';
                case {'biofilm'}
                    b = model.getComponentNames();
                    index = find(strcmpi(b, 'biofilm')) - 3;
                    fn = 'b';
                case {'calcite'}
                    c = model.getComponentNames();
                    index = find(strcmpi(c, 'calcite')) - 4;
                    fn = 'c';
                otherwise
                    [fn, index] = ...
                     getVariableField@TwoPhaseOilWaterModel(model, ...
                                                        name, varargin{:});
            end
        end

        function names = getComponentNames(model)
            names = getComponentNames@TwoPhaseOilWaterModel(model);
            names{end+1} = 'microorganism';
            names{end+1} = 'oxygen';
            names{end+1} = 'urea';
            names{end+1} = 'biofilm';
            names{end+1} = 'calcite';
        end

        function scaling = getScalingFactorsCPR(model, problem, names, ...
                                                                    solver)
            nNames = numel(names);
            scaling = cell(nNames, 1);
            handled = false(nNames, 1);
            for iter = 1 : nNames
                name = lower(names{iter});
                switch name
                    case 'microorganism'
                        s = 0;
                    case 'oxygen'
                        s = 0;
                    case 'urea'
                        s = 0; 
                    case 'biofilm'
                        s = 0;
                    case 'calcite'
                        s = 0;
                    otherwise
                        continue
                end
                sub = strcmpi(problem.equationNames, name);
                scaling{iter} = s;
                handled(sub) = true;
            end
            if ~all(handled)
                other = getScalingFactorsCPR@ThreePhaseBlackOilModel( ...
                                  model, problem, names(~handled), solver);
                [scaling{~handled}] = other{:};
            end
        end

        function [eq, src] = addComponentContributions(model, cname, ...
                                                 eq, component, src, force)
            if isempty(force)
                return
            end
            c = model.getProp(force, cname);
            cells = src.sourceCells;
            switch lower(cname)
              case {'microorganism'}
                qW = src.phaseMass{1} ./ model.fluid.rhoWS;
                isInj = qW > 0;
                qC = (isInj .* c + ~isInj .* component(cells)) .* qW;  
              case {'oxygen'}
                qW = src.phaseMass{1} ./ model.fluid.rhoWS;
                isInj = qW > 0;
                qC = (isInj .* c + ~isInj .* component(cells)) .* qW;
              case {'urea'}
                qW = src.phaseMass{1} ./ model.fluid.rhoWS;
                isInj = qW > 0;
                qC = (isInj .* c + ~isInj .* component(cells)) .* qW;
              case {'biofilm'}
                qC = 0;
              case {'calcite'}
                qC = 0;
              otherwise
                error(['Unknown component ''', cname, ...
                                               '''. BC not implemented.']);
            end
            eq(cells) = eq(cells) - qC;
            src.components{end+1} = qC;
        end

        function [compEqs, compSrc, eqNames, wellSol] = ...
                  getExtraWellContributions(model, well, wellSol0, ...
                      wellSol, q_s, bh, packed, qMass, qVol, dt, iteration)
            [compEqs, compSrc, eqNames, wellSol] = ...
             getExtraWellContributions@TwoPhaseOilWaterModel(model, ...
             well, wellSol0, wellSol, q_s, bh, packed, qMass, qVol, ...
                                                            dt, iteration);
            f = model.fluid;
            wix = 1; 
            cqWs = qMass{wix} ./ f.rhoWS; 
            q_sw = q_s{wix};  
            isInj = (cqWs > 0);
            qWout = sum(cqWs(isInj));
            if (q_sw < 0) 
                qWout = qWout - q_sw;
            end
            conccompCtrl = model.getProp(well.W, 'microorganism');
            pix = strcmpi(model.getComponentNames(), 'microorganism');
            conccompRes = packed.components{pix};
            qC = sum(-conccompRes(~isInj) .* cqWs(~isInj));
            if (q_sw > 0) 
                qC = qC + conccompCtrl * q_sw;
            end               
            cqC   = conccompRes .* cqWs;
            if any(isInj)
                cqC(isInj) = (qC ./ qWout) .* cqWs(isInj);
            end                               
            compSrc{end+1} = cqC;
            conccompCtrl = model.getProp(well.W, 'oxygen');
            pix = strcmpi(model.getComponentNames(), 'oxygen');
            conccompRes = packed.components{pix};
            qC = sum(-conccompRes(~isInj) .* cqWs(~isInj));
            if (q_sw > 0) 
                qC = qC + conccompCtrl * q_sw;
            end               
            cqC   = conccompRes .* cqWs;
            if any(isInj)
                cqC(isInj) = (qC ./ qWout) .* cqWs(isInj);
            end                               
            compSrc{end+1} = cqC;
            conccompCtrl = model.getProp(well.W, 'urea');
            pix = strcmpi(model.getComponentNames(), 'urea');
            conccompRes = packed.components{pix};
            qC = sum(-conccompRes(~isInj) .* cqWs(~isInj));
            if (q_sw > 0) 
                qC = qC + conccompCtrl * q_sw;
            end               
            cqC   = conccompRes .* cqWs;
            if any(isInj)
                cqC(isInj) = (qC ./ qWout) .* cqWs(isInj);
            end                               
            compSrc{end+1} = cqC;
            pix = strcmpi(model.getComponentNames(), 'biofilm');
            conccompRes = packed.components{pix};
            compSrc{end+1} = conccompRes .* 0;
            pix = strcmpi(model.getComponentNames(), 'calcite');
            conccompRes = packed.components{pix};
            compSrc{end+1} = conccompRes .* 0;
        end
    end
end