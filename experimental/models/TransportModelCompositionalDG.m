classdef TransportModelCompositionalDG < TransportModelDG
   
    properties
    end
    
    methods
        function model = TransportModelCompositionalDG(parent, varargin)
            model = model@TransportModelDG(parent, varargin{:});
        end
        
        function [state, report] = updateState(model, state, problem, dx, drivingForces)
            
            parent = model.parentModel;
            state0 = state;
            
            var0 = problem.primaryVariables;
            vars = var0;
            removed = false(size(vars));
            
            wix = strcmpi(vars, 'sW');
            if any(wix)
                state = model.updateStateFromIncrement(state, dx{wix}, problem, 'sW', inf, model.dsMaxAbs);
                removed(wix) = true;
                vars = vars(~wix);
            end
            
            % Components
            cnames = parent.EOSModel.fluid.names;
            ncomp = numel(cnames);
            ok = false(ncomp, 1);

            z    = model.getProp(state, 'components');
            zdof = model.getProp(state, 'componentsdof');
            rm = 0;

            for i = 1:ncomp
                name = [lower(cnames{i}), 'dof'];
                cix = strcmpi(var0, name);
                if any(cix)
                    zdof0 = zdof(:, i);

                    dz = dx{cix};
                    if isfinite(parent.dzMaxAbs)
                        dz = sign(dz).*min(abs(dz), parent.dzMaxAbs);
                    end
                    zdof(:, i) = min(max(zdof0 + dz, 0), 1);

                    ok(i) = true;
                    [vars, ix] = model.stripVars(vars, {name});
                    removed(~removed) = removed(~removed) | ix;

                    rm = rm - (zdof(:, i) - zdof0);
                    
                end
            end
            
            if any(ok)
                % We had components as active variables somehow
                assert(nnz(~ok) == 1)
                z(:, ~ok) = min(max(z(:, ~ok) + rm, 0), 1);
                z = bsxfun(@rdivide, z, sum(z, 2));
                state.components = z;
                if model.water
                    v  = 1 - state.s(:, 1);
                    v0 = 1 - state0.s(:, 1);
                else
                    [v, v0] = deal(1);
                end

                state.dz = computeChange(z, state0.components, v, v0);
            end

            % Parent class handles almost everything for us
            problem.primaryVariables = vars;
            dx(removed) = [];
            [state, report] = updateState@ThreePhaseCompositionalModel(model, state, problem, dx, drivingForces);
            
            if problem.iterationNo == 1
                state.switched = false(model.G.cells.num, 1);
                state.switchCount = zeros(model.G.cells.num, 1);
            end
            twoPhase0 = state.L < 1 & state.L > 0;
            % Update saturations etc using flash routine
            state = model.computeFlash(state, problem.dt, problem.iterationNo);
            twoPhase = state.L < 1 & state.L > 0;
            switched = twoPhase0 ~= twoPhase;
            
            dispif(model.verbose > 1, '%d gas, %d oil, %d two-phase\n', nnz(state.L == 0), nnz(state.L == 1), nnz(twoPhase));
            state.switchCount = state.switchCount + double(switched);
            
            minz = model.EOSModel.minimumComposition;
            state.components = ensureMinimumFraction(state.components, minz);
            state.x = ensureMinimumFraction(state.x, minz);
            state.y = ensureMinimumFraction(state.y, minz);
        end
        
         %-----------------------------------------------------------------%
        function [state, names, origin] = getStateAD(model, state, init)
            if nargin < 3
                init = true;
            end
            parent = model.parentModel;
            % Get the AD state for this model
            [basevars, dofbasenames, basenames, baseorigin] = model.getPrimaryVariables(state);
            % Find saturations
            isS = false(size(basevars));
            nph = parent.getNumberOfPhases();
            phase_variable_index = zeros(nph, 1);
            for i = 1:numel(basevars)
                [f, ix] = model.getVariableField(dofbasenames{i}, false);
                if strcmp(f, 'sdof') || strcmpi(dofbasenames{i}, 'xdof')
                    isS(i) = true;
                    phase_variable_index(ix) = i;
                end
            end
            % Figure out saturation logic
            isP    = strcmp(dofbasenames, 'pressuredof');
            vars   = basevars;
            names  = dofbasenames;
            origin = baseorigin;
            useTotalSaturation = strcmpi(model.formulation, 'totalSaturation') ...
                                    && sum(isS) == nph - 1;
            useTotalSaturation  = useTotalSaturation || strcmpi(class(parent), 'GenericOverallCompositionModel');
%             assert(useTotalSaturation, 'DG currently only supports total saturation formulation!');
            if useTotalSaturation
                % Replace pressure with total saturation
                replacement = 'sTdof';
                sTdof = model.getProp(state, replacement);
                % Replacing
                vars{isP} = sTdof;
                names{isP} = replacement;
                origin{isP} = class(model);
            else
                % Remove pressure and skip saturation closure
                vars = vars(~isP);
                names = names(~isP);
                origin = origin(~isP);
            end
            if init
                [vars{:}] = model.AutoDiffBackend.initVariablesAD(vars{:});
            end
            if useTotalSaturation
                basevars(~isP) = vars(~isP);
            else
                basevars(~isP) = vars;
            end
            % Let parent model handle state initialization
            state = model.initStateAD(state, basevars, dofbasenames, baseorigin);
            fixedSat = false;
            for i = 1:numel(dofbasenames)
                if any(strcmpi(basenames{i}, {'sw', 'so', 'sg', 'x'}))
                    if fixedSat
                        continue
                    end
                    basenames{i} = 's';
                    dofbasenames{i}  = 'sdof';
                    fixedSat = true;
                end
                v     = model.getProp(state, dofbasenames{i});
                vm    = model.discretization.getCellMean(state, value(v));
                state = model.setProp(state, basenames{i}, vm); 
            end
            if useTotalSaturation
                % Set total saturation as well
                sTdof       = vars{isP};
                state.sTdof = sTdof;
                % Evaluate at cell cubature points
                cellValue         = model.discretization.evaluateProp(state, sTdof, 'cell');
                state.cellStateDG = model.setProp(state.cellStateDG, 'sT', cellValue);
                % Evaluate mean
                cellMean          = model.discretization.getCellMean(state, sTdof);
                state.wellStateDG = model.setProp(state.wellStateDG, 'sT', cellMean);
                % Evaluate at face cubature points
                faceValue         = model.discretization.evaluateProp(state, sTdof, 'face');
                state.faceStateDG = model.setProp(state.faceStateDG, 'sT', faceValue);
                % Set mean in state
                state = model.setProp(state, 'sT', cellMean);
            end
        end
        
        function state = initStateAD(model, state, vars, names, origin)
            
            pmodel = model.parentModel;
            
            isP   = strcmp(names, 'pressuredof');
            state = model.setProp(state, 'pressuredof', vars{isP});
            p     = model.discretization.getCellMean(state, vars{isP});
            state = model.setProp(state, 'pressure', p);
            
            pCell = model.discretization.evaluateProp(state, vars{isP}, 'cell');
            pFace = model.discretization.evaluateProp(state, vars{isP}, 'face');
            
            removed = isP;
            
            cnames = pmodel.EOSModel.fluid.names;
            ncomp = numel(cnames);
            zdof = cell(1, ncomp);
            z_end = model.discretization.getFillSat(state);
            for i = 1:ncomp
                name = cnames{i};
                sub = strcmp(names, [name, 'dof']);
                if any(sub)
                    zdof{i} = vars{sub};
                    z_end = z_end - zdof{i};
                    removed(sub) = true;
                else
                    fill = i;
                end
            end
            zdof{fill} = z_end;
            state = model.setProp(state, 'componentsdof', zdof);
            z     = model.discretization.getCellMean(state, zdof);
            state = model.setProp(state, 'components', z);
            
            zCell = model.discretization.evaluateProp(state, zdof, 'cell');
            zFace = model.discretization.evaluateProp(state, zdof, 'face');
            
            TCell = model.discretization.evaluateProp(state, state.Tdof, 'cell');
            
            if isa(state.pressure, 'ADI') || isa(zdof{1}, 'ADI')
                [state.x, state.y, state.L] = pmodel.EOSModel.getPhaseFractionAsADI(state, state.pressure, state.T, state.components);
                [xCell, yCell, LCell] = pmodel.EOSModel.getPhaseFractionAsADI(state, pCell, TCell, zCell);
            end
            if ~isempty(pmodel.FacilityModel)
                % Select facility model variables and pass them off to attached
                % class.
                fm = class(pmodel.FacilityModel);
                isF = strcmp(origin, fm);
                state = pmodel.FacilityModel.initStateAD(state, vars(isF), names(isF), origin(isF));
                removed = removed | isF;
            end
            if pmodel.water
                isWater = strcmp(names, 'water');
                sW = vars{isWater};
                removed(isWater) = true;
                offset = 1;
            else
                offset = 0;
            end
            % Set up state with remaining variables
            state = pmodel.initStateFunctionContainers(state);
            %             state = initStateAD@GenericOverallCompositionModel(pmodel, state, vars(~removed), names(~removed), origin(~removed));
            
            % Now that props have been set up, we can get the saturations
            Z = pmodel.getProps(state, 'PhaseCompressibilityFactors');
            Z_L = Z{offset+1};
            Z_V = Z{offset+2};
            
            L = state.L;
            volL = L.*Z_L;
            volV = (1-L).*Z_V;
            volT = volL + volV;
            sL = volL./volT;
            sV = volV./volT;
            
            if pmodel.water
                [pureLiquid, pureVapor] = model.getFlag(state);
                void = 1 - sW;
                sL = sL.*void;
                sV = sV.*void;
                [sL, sV] = model.setMinimumTwoPhaseSaturations(state, sW, sL, sV, pureVapor, pureLiquid);

                s = {sW, sL, sV};
            else
                s = {sL, sV};
            end
            state = model.setProp(state, 's', s);
            state = model.evaluateBaseVariables(state);
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
