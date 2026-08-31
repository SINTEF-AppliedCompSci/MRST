classdef MICPModel < TwoPhaseOilWaterModel
% Script to implement the MICP model.
%
% This script is modified from a file in the MATLAB Reservoir Simulation
% Toolbox (MRST), see
%   mrst/modules/ad-eor/models/OilWaterPolymerModel.m
%
% We refer to that script for a complete commented version of the file.

%{
Partial copyright 2009-2021, SINTEF Digital, Mathematics & Cybernetics.
Partial copyright 2021-2026, NORCE Research AS, Computational Geosciences and
Modeling.

This file is part of the ad-micp module.

ad-micp is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

ad-micp is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this file. If not, see <http://www.gnu.org/licenses/>.
%}
    properties
        % Initial porosity before biofilm and calcite accumulation
        referencePorosity

        % Minimum porosity used to prevent numerical singularities
        minimumPorosity = 1e-8
    end

    methods
        function model = MICPModel(G, rock, fluid, varargin)
            model = model@TwoPhaseOilWaterModel(G, rock, fluid);
            model.referencePorosity = rock.poro;
            model = merge_options(model, varargin{:});

            assert(isscalar(model.minimumPorosity) && ...
                model.minimumPorosity > 0, ...
                'minimumPorosity must be a positive scalar.');

            assert( ...
                isscalar(model.referencePorosity) || ...
                numel(model.referencePorosity) == model.G.cells.num, ...
                ['referencePorosity must be scalar or contain one ' ...
                 'value per active grid cell.']);

            assert(all(model.referencePorosity(:) > ...
                model.minimumPorosity), ...
                ['referencePorosity must be larger than ' ...
                 'minimumPorosity in every cell.']);

            assert(isfield(model.fluid, 'kmin') && ...
                all(model.fluid.kmin(:) > 0), ...
                'fluid.kmin must be defined and positive.');
        end

        function [problem, state] = getEquations( ...
                model, state0, state, dt, drivingForces, varargin)

            [problem, state] = equationsMICP( ...
                state0, state, model, dt, drivingForces, varargin{:});
        end

        function state = validateState(model, state)
            state = validateState@TwoPhaseOilWaterModel(model, state);

            % Components must be present
            model.checkProperty( ...
                state, 'Microorganism', model.G.cells.num, 1);
            model.checkProperty( ...
                state, 'Oxygen', model.G.cells.num, 1);
            model.checkProperty( ...
                state, 'Urea', model.G.cells.num, 1);
            model.checkProperty( ...
                state, 'Biofilm', model.G.cells.num, 1);
            model.checkProperty( ...
                state, 'Calcite', model.G.cells.num, 1);
        end

        function [state, report] = updateState( ...
                model, state, problem, dx, drivingForces)

            % Store the variables from the previous iteration temporarily
            % for use in the convergence criteria
            microorganismPrevious = ...
                model.getProp(state, 'microorganism');
            oxygenPrevious = model.getProp(state, 'oxygen');
            ureaPrevious = model.getProp(state, 'urea');
            biofilmPrevious = model.getProp(state, 'biofilm');
            calcitePrevious = model.getProp(state, 'calcite');

            [state, report] = ...
                updateState@TwoPhaseOilWaterModel( ...
                model, state, problem, dx, drivingForces);

            % Limit the dissolved variables to [0, cmax]. For the
            % microorganisms, we set this value equal to the maximum
            % biomass density found in the literature (105 kg/m^3). For
            % oxygen and urea, we set this value to the maximum injected
            % concentration. The combined biofilm and calcite volume is
            % limited by the reference porosity minus the minimum porosity.
            state = model.limitComponent( ...
                state, 'microorganism', 'm_prev', ...
                microorganismPrevious, 105);

            state = model.limitComponent( ...
                state, 'oxygen', 'o_prev', ...
                oxygenPrevious, model.fluid.Comax);

            state = model.limitComponent( ...
                state, 'urea', 'u_prev', ...
                ureaPrevious, model.fluid.Cumax);

            state = model.limitSolidComponents( ...
                state, biofilmPrevious, calcitePrevious);
        end

        function [state, report] = updateAfterConvergence( ...
                model, state0, state, dt, drivingForces)

            [state, report] = ...
                updateAfterConvergence@TwoPhaseOilWaterModel( ...
                model, state0, state, dt, drivingForces);

            % Remove the temporary fields used for convergence
            temporaryFields = { ...
                'm_prev', ...
                'o_prev', ...
                'u_prev', ...
                'b_prev', ...
                'c_prev'};

            fieldExists = false(size(temporaryFields));

            for fieldIndex = 1 : numel(temporaryFields)
                fieldExists(fieldIndex) = ...
                    isfield(state, temporaryFields{fieldIndex});
            end

            temporaryFields = temporaryFields(fieldExists);

            if ~isempty(temporaryFields)
                state = rmfield(state, temporaryFields);
            end
        end

        function names = getComponentNames(model)
            names = getComponentNames@TwoPhaseOilWaterModel(model);
            names{end + 1} = 'microorganism';
            names{end + 1} = 'oxygen';
            names{end + 1} = 'urea';
            names{end + 1} = 'biofilm';
            names{end + 1} = 'calcite';
        end

        function scaling = getScalingFactorsCPR( ...
                model, problem, names, solver)

            numberOfNames = numel(names);
            scaling = cell(numberOfNames, 1);
            handled = false(numberOfNames, 1);

            micpComponents = { ...
                'microorganism', ...
                'oxygen', ...
                'urea', ...
                'biofilm', ...
                'calcite'};

            for nameIndex = 1 : numberOfNames
                name = names{nameIndex};

                if any(strcmpi(name, micpComponents))
                    scaling{nameIndex} = 0;
                    handled(nameIndex) = true;
                end
            end

            if ~all(handled)
                % Get rest of scaling factors
                other = getScalingFactorsCPR@ThreePhaseBlackOilModel(model, problem, names(~handled), solver);
                [scaling{~handled}] = other{:};
            end
        end

        function [equation, source] = addComponentContributions( ...
                model, componentName, equation, component, source, force)

            if isempty(force)
                return
            end

            controlledConcentration = ...
                model.getProp(force, componentName);
            cells = source.sourceCells;

            switch lower(componentName)
                case {'microorganism', 'oxygen', 'urea'}
                    waterRate = ...
                        source.phaseMass{1} ./ model.fluid.rhoWS;
                    isInjection = waterRate > 0;

                    componentRate = ( ...
                        isInjection .* controlledConcentration + ...
                        ~isInjection .* component(cells)) .* waterRate;

                case {'biofilm', 'calcite'}
                    componentRate = 0;

                otherwise
                    error( ...
                        'Unknown component ''%s''. BC not implemented.', ...
                        componentName);
            end

            equation(cells) = equation(cells) - componentRate;
            source.components{end + 1} = componentRate;
        end

        function [fieldName, index] = getVariableField( ...
                model, name, varargin)

            % Get the index/name mapping for the model
            switch lower(name)
                case 'microorganism'
                    fieldName = 'm';
                    index = 1;

                case 'oxygen'
                    fieldName = 'o';
                    index = 1;

                case 'urea'
                    fieldName = 'u';
                    index = 1;

                case 'biofilm'
                    fieldName = 'b';
                    index = 1;

                case 'calcite'
                    fieldName = 'c';
                    index = 1;

                otherwise
                    [fieldName, index] = ...
                        getVariableField@TwoPhaseOilWaterModel( ...
                        model, name, varargin{:});
            end
        end

        function [componentEquations, componentSources, ...
                equationNames, wellSol] = getExtraWellContributions( ...
                model, well, wellSol0, wellSol, q_s, bh, packed, ...
                qMass, qVol, dt, iteration)

            [componentEquations, componentSources, ...
                equationNames, wellSol] = ...
                getExtraWellContributions@TwoPhaseOilWaterModel( ...
                model, well, wellSol0, wellSol, q_s, bh, packed, ...
                qMass, qVol, dt, iteration);

            fluid = model.fluid;
            componentNames = model.getComponentNames();

            transportedComponents = { ...
                'microorganism', ...
                'oxygen', ...
                'urea'};

            immobileComponents = { ...
                'biofilm', ...
                'calcite'};

            waterPhaseIndex = 1;
            completionWaterRates = ...
                qMass{waterPhaseIndex} ./ fluid.rhoWS;
            surfaceWaterRate = q_s{waterPhaseIndex};
            isInjection = completionWaterRates > 0;

            totalOutgoingWaterRate = ...
                sum(completionWaterRates(isInjection));

            if surfaceWaterRate < 0
                totalOutgoingWaterRate = ...
                    totalOutgoingWaterRate - surfaceWaterRate;
            end

            if any(isInjection)
                outgoingWaterRateValue = ...
                    value(totalOutgoingWaterRate);

                assert(all(outgoingWaterRateValue(:) > 0), ...
                    ['Injection completions require a positive ' ...
                     'total outgoing water rate.']);
            end

            for componentIndex = 1 : numel(transportedComponents)
                componentName = ...
                    transportedComponents{componentIndex};

                componentMatches = ...
                    strcmpi(componentNames, componentName);

                assert(nnz(componentMatches) == 1, ...
                    sprintf( ...
                    ['Expected exactly one component named ''%s'', ' ...
                     'but found %d.'], ...
                    componentName, nnz(componentMatches)));

                packedComponentIndex = find( ...
                    componentMatches, 1);

                controlledConcentration = ...
                    model.getProp(well.W, componentName);

                reservoirConcentration = ...
                    packed.components{packedComponentIndex};

                totalComponentRate = sum( ...
                    -reservoirConcentration(~isInjection) .* ...
                    completionWaterRates(~isInjection));

                if surfaceWaterRate > 0
                    totalComponentRate = totalComponentRate + ...
                        controlledConcentration * surfaceWaterRate;
                end

                completionComponentRates = ...
                    reservoirConcentration .* completionWaterRates;

                if any(isInjection)
                    completionComponentRates(isInjection) = ...
                        (totalComponentRate ./ ...
                        totalOutgoingWaterRate) .* ...
                        completionWaterRates(isInjection);
                end

                componentSources{end + 1} = ...
                    completionComponentRates;
            end

            for componentIndex = 1 : numel(immobileComponents)
                componentName = ...
                    immobileComponents{componentIndex};

                componentMatches = ...
                    strcmpi(componentNames, componentName);

                assert(nnz(componentMatches) == 1, ...
                    sprintf( ...
                    ['Expected exactly one component named ''%s'', ' ...
                     'but found %d.'], ...
                    componentName, nnz(componentMatches)));

                packedComponentIndex = find( ...
                    componentMatches, 1);

                reservoirConcentration = ...
                    packed.components{packedComponentIndex};

                componentSources{end + 1} = ...
                    reservoirConcentration .* 0;
            end
        end
    end

    methods (Access = private)
        function state = limitComponent( ...
                model, state, componentName, previousFieldName, ...
                previousValue, upperBound)

            component = model.getProp(state, componentName);
            component = min(component, upperBound);
            component = max(component, 0);

            state = model.setProp( ...
                state, componentName, component);
            state.(previousFieldName) = previousValue;
        end

        function state = limitSolidComponents( ...
                model, state, biofilmPrevious, calcitePrevious)

            biofilm = model.getProp(state, 'biofilm');
            calcite = model.getProp(state, 'calcite');

            biofilm = max(biofilm, 0);
            calcite = max(calcite, 0);

            maximumSolidVolume = max( ...
                model.referencePorosity - model.minimumPorosity, 0);

            if isscalar(maximumSolidVolume)
                maximumSolidVolume = ...
                    repmat(maximumSolidVolume, size(biofilm));
            else
                maximumSolidVolume = ...
                    reshape(maximumSolidVolume, size(biofilm));
            end

            totalSolidVolume = biofilm + calcite;
            constrainedCells = ...
                totalSolidVolume > maximumSolidVolume;

            if any(constrainedCells)
                scalingFactor = ...
                    maximumSolidVolume(constrainedCells) ./ ...
                    totalSolidVolume(constrainedCells);

                biofilm(constrainedCells) = ...
                    biofilm(constrainedCells) .* scalingFactor;

                calcite(constrainedCells) = ...
                    calcite(constrainedCells) .* scalingFactor;
            end

            state = model.setProp(state, 'biofilm', biofilm);
            state = model.setProp(state, 'calcite', calcite);
            state.b_prev = biofilmPrevious;
            state.c_prev = calcitePrevious;
        end
    end
end
