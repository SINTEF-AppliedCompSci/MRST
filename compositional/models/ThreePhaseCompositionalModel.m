classdef ThreePhaseCompositionalModel < ReservoirModel
    % Base class for compositional models
    %
    % SYNOPSIS:
    %   model = ThreePhaseCompositionalModel(G, rock, fluid, compFluid)
    %
    % DESCRIPTION:
    %   This is the base class for several compositional models in MRST. It
    %   contains common functionality and is not intended for direct use.
    %
    % PARAMETERS:
    %   G         - Grid structure
    %   rock      - Rock structure for the reservoir
    %   fluid     - The flow fluid, containing relative permeabilities,
    %               surface densities and flow properties for the
    %               aqueous/water phase (if present)
    %   compFluid - CompositionalFluid instance describing the species
    %               present.
    %
    % RETURNS:
    %  model - Initialized class instance
    %
    % SEE ALSO:
    %   `NaturalVariablesCompositionalModel`, `OverallCompositionCompositionalModel`

    properties
        EOSModel % EquationOfState model used for PVT and phase behavior
        EOSNonLinearSolver % NonLinearSolver used to solve EOS problems that appear
        dzMaxAbs % Maximum allowable change in any mole fraction
        incTolPressure % Relative increment tolerance for pressure
        incTolComposition % Increment tolerance for composition
        useIncTolComposition % If true, use increment tolerance for composition. Otherwise, use mass-balance.
        fugacityTolerance % Tolerance for fugacity equality (in units 1/barsa)
    end
    
    methods
        function model = ThreePhaseCompositionalModel(G, rock, fluid, compFluid, varargin)
            model = model@ReservoirModel(G, rock, fluid);
            
            if isa(compFluid, 'CompositionalFluid')
                model.EOSModel = EquationOfStateModel(G, compFluid);
            elseif isa(compFluid, 'EquationOfStateModel')
                model.EOSModel = compFluid;
            end
            
            model.EOSModel.verbose = false;
            model.EOSNonLinearSolver = getDefaultFlashNonLinearSolver();
            
            model.nonlinearTolerance = 0.01;
            model.incTolPressure = 1e-3;
            model.useIncTolComposition = false;
            model.incTolComposition = 1e-3;
            
            model.fugacityTolerance = 1e-3;
            
            model.water = true;
            model.oil = true;
            model.gas = true;
            
            model.dzMaxAbs = 0.1;
            model.dsMaxAbs = 0.1;
            
            model.minimumPressure = 0;
            model.dpMaxRel = 0.25;
            model = merge_options(model, varargin{:});
        end
        
        function [problem, state] = getEquations(model, varargin)
            error('Superclass not intended for direct use!');
        end
        
        function names = getComponentNames(model)
            names = getComponentNames@ReservoirModel(model);
            names = horzcat(names, model.EOSModel.fluid.names);
        end
        
        function scaling = getComponentScaling(model, state)
            wL = state.s(:, 1+model.water);
            wV = state.s(:, 2+model.water);
            wT = wL + wV;

            wL = wL./wT;
            wV = wV./wT;
            
            if isfield(state, 'rho')
                rhoO = state.rho(:, 1+model.water);
                rhoG = state.rho(:, 2+model.water);
            else
                rhoO = model.fluid.rhoOS;
                rhoG = model.fluid.rhoGS;
            end
            wL(wT == 0) = state.L(wT == 0);
            wV(wT == 0) = 1 - state.L(wT == 0);

            scaling = wL.*double(rhoO) + wV.*double(rhoG);
            scaling(wT  < 1e-3) = 1e6;
        end

        function [fn, index] = getVariableField(model, name)
            switch(lower(name))
                case {'z', 'components'}
                    % Overall mole fraction
                    fn = 'components';
                    index = ':';
                case {'x', 'liquidmf'}
                    % Liquid phase mole fraction
                    fn = 'x';
                    index = ':';
                case {'y', 'vapormf'}
                    % Vapor phase mole fraction
                    fn = 'y';
                    index = ':';
                case {'l', 'liquid'}
                    % Liquid mole fraction
                    fn = 'L';
                    index = 1;
                case {'z_l'}
                    % Liquid compressibility
                    fn = 'Z_L';
                    index = 1;
                case {'z_v'}
                    % Vapor compressibility
                    fn = 'Z_V';
                    index = 1;
                otherwise
                    sub = strcmpi(model.getComponentNames(), name);
                    if any(sub)
                        fn = 'components';
                        index = find(sub);
                    else
                        % This will throw an error for us
                        [fn, index] = getVariableField@ReservoirModel(model, name);
                    end
            end
        end

        function state = validateState(model, state)
            state = validateState@ReservoirModel(model, state);
            ncell = model.G.cells.num;
            ncomp = model.EOSModel.fluid.getNumberOfComponents();
            model.checkProperty(state, 'Components', [ncell, ncomp], [1, 2]);
            assert(all(max(model.getProp(state, 'Components')) <= 1), ...
                'Initial mole fractions are larger than unity.')
            T = model.getProp(state, 'Temperature');
            if numel(T) == 1
                % Expand single temperature to all grid cells
                fn = model.getVariableField('Temperature');
                state = rmfield(state, fn);
                state = model.setProp(state, 'Temperature', repmat(T, ncell, 1));
            end
            model.checkProperty(state, 'Temperature', [ncell, 1], [1, 2]);
            if ~isfield(state, 'x') || ~isfield(state, 'K')
                state.components = ensureMinimumFraction(state.components, model.EOSModel.minimumComposition);
                state = model.computeFlash(state, inf);
            end
            if isfield(state, 'wellSol') && ~isfield(state.wellSol, 'components')
                for i = 1:numel(state.wellSol)
                    state.wellSol(i).components = [];
                end
            end
        end

        function [state, report] = updateState(model, state, problem, dx, drivingForces)
            p0 = state.pressure;
            [state, report] = updateState@ReservoirModel(model, state, problem, dx, drivingForces);
            range = max(p0) - min(p0);
            if range == 0
                range = 1;
            end
            state.dpRel = (state.pressure - p0)./range;
            state.dpAbs = state.pressure - p0;
        end
        
        function state = computeFlash(model, state, dt, iteration)
            % Flash a state with the model's EOS.
            if nargin < 4
                iteration = 1;
            end
            state0 = state;
            if iteration == 1
                state.eos.iterations = 0;
                state.eos.cellcount = 0;
                state = model.EOSModel.validateState(state);
            end
            [state, report] = model.EOSNonLinearSolver.solveTimestep(state, dt, model.EOSModel);

            if ~isempty(report)
                state.eos.iterations = state.eos.iterations + report.StepReports{1}.Iterations;

                if ~report.StepReports{1}.Converged
                    state = model.EOSModel.updateAfterConvergence(state0, state, dt, struct());
                    disp   ('********************************');
                    disp   ('*    Flash did not converge    *');
                    disp   ('********************************');
                    fprintf('* Final residuals after %d iterations:\n', report.StepReports{1}.Iterations);
                    for i = 1:numel(report.StepReports{1}.NonlinearReport{end}.Residuals)
                        fprintf('* %s: %1.4g \n', model.EOSModel.fluid.names{i},....
                                report.StepReports{1}.NonlinearReport{end}.Residuals(i))
                    end
                end
            end
            L = state.L;
            Z_L = state.Z_L;
            Z_V = state.Z_V;

            if model.water
                sW = state.s(:, 1);
                sL = (1 - sW).*L.*Z_L./(L.*Z_L + (1-L).*Z_V);
                sV = 1 - sW - sL;
                state.s = [sW, sL, sV];
            else
                sL = L.*Z_L./(L.*Z_L + (1-L).*Z_V);
                sV = 1 - sL;
                state.s = [sL, sV];
            end
        end
    
        
        function [v_eqs, tolerances, names] = getConvergenceValues(model, problem, varargin)
            [v_wells, tol_wells, names_wells, is_well] = ...
                model.FacilityModel.getFacilityConvergenceValues(problem);

            v_eqs = norm(problem, inf);
            % Check components
            [v_comp, tol_comp, names_comp, is_comp] = model.getComponentConvergenceValues(problem);
            % Check fugacity
            [v_f, tol_f, names_f, is_f] = model.getFugacityConvergenceValues(problem);
            % Remaining values use nonlinear tolerances
            rest = ~(is_well | is_f | is_comp);
            v_rest = v_eqs(rest);
            tol_rest = repmat(model.nonlinearTolerance, size(v_rest));
            names_rest = problem.equationNames(rest);
            % Define tolerances, values and names
            tolerances = [tol_comp, tol_rest, tol_f, tol_wells];
            v_eqs = [v_comp, v_rest, v_f, v_wells];
            names = horzcat(names_comp, names_rest, names_f, names_wells);
            % Pressure tolerances
            if any(strcmpi(problem.primaryVariables, 'pressure'))
                if isfield(problem.state, 'dpRel') && problem.iterationNo > 1
                    dp = norm(problem.state.dpRel, inf);
                else
                    dp = inf;
                end
                if isfinite(model.incTolPressure)
                    isp = strcmpi(names, 'pressure');
                    if any(isp)
                        v_eqs(isp) = dp;
                        tolerances(isp) = model.incTolPressure;
                        names{isp} = 'dPressure';
                    else
                        v_eqs = [dp, v_eqs];
                        tolerances = [model.incTolPressure, tolerances];
                        names = ['deltaP', names];
                    end
                elseif isa(model, 'PressureNaturalVariablesModel')
                    pRes = norm(problem.b, inf);
                    if isempty(pRes)
                        pRes = inf;
                    end
                    v_eqs = [pRes, v_eqs];
                    tolerances = [model.nonlinearTolerance, tolerances];
                    names = ['Pressure', names];
                end
            end
        end
        
        function state = storeDensities(model, state, rhoW, rhoO, rhoG)
            % Densities
            isActive = model.getActivePhases();

            state.rho = zeros(model.G.cells.num, sum(isActive));
            rho = {double(rhoW), double(rhoO), double(rhoG)};
            state = model.setPhaseData(state, rho, 'rho');
        end
        
        function [compEqs, compSrc, eqNames, wellSol] = getExtraWellContributions(model, well, wellSol0, wellSol, q_s, bh, packed, qMass, qVol, dt, iteration)
            [compEqs, compSrc, eqNames, wellSol] = getExtraWellContributions@ReservoirModel(model, well, wellSol0, wellSol, q_s, bh, packed, qMass, qVol, dt, iteration);
            names = model.getComponentNames();
            ncomp = numel(names);
            
            
            z_well = model.getProp(well.W, 'components');
            mf_well = model.EOSModel.getMassFraction(z_well);
            
            lix = model.water + 1;
            vix = lix + 1;
            cqLs = qMass{lix};
            cqVs = qMass{vix};
            ncell = numel(double(cqLs));

            N = numel(compSrc);
            compSrc = [compSrc, cell(1, ncomp)];
            
            wellSol.components = zeros(ncell, ncomp);
            for cNo = 1:ncomp
                Z_well = mf_well(cNo);
                X_res = packed.components{cNo}{1};
                Y_res = packed.components{cNo}{2};

                injO = cqLs > 0;
                injG = cqVs > 0;
                % Account for both phases.
                q_i = (cqLs.*injO + cqVs.*injG).*Z_well ...
                       + ~injO.*X_res.*cqLs + ~injG.*Y_res.*cqVs;


                compSrc{N+cNo} = q_i;
                wellSol.components(:, cNo) = double(q_i);
            end
        end

        function [eq, src] = addComponentContributions(model, cname, eq, component, src, force)
            if isempty(force)
                return
            end
            cnames = model.getComponentNames();
            sub = strcmpi(cnames, cname);
            if any(sub)
                cells = src.sourceCells;
                if isfield(force, 'componentMass')
                    qC = force.componentMass(:, sub);
                else
                    if isfield(force, 'x')
                        assert(isfield(force, 'y'))
                        [x_bc, y_bc] = model.getProps(force, 'x', 'y');
                        massFractions = {model.EOSModel.getMassFraction(x_bc), ...
                                         model.EOSModel.getMassFraction(y_bc)};
                    else
                        z_bc = model.getProp(force, 'components');
                        mf_bc = model.EOSModel.getMassFraction(z_bc);
                        
                        massFractions = {mf_bc, mf_bc};
                    end
                    qC = zeros(size(cells));
                    for ph = 1:2
                        ix = ph + model.water;
                        q_ph = src.phaseMass{ix};
                        inj = q_ph > 0;

                        qC = qC + ~inj.*component{ph}(cells).*q_ph ...
                                +  inj.*massFractions{ph}(:, sub).*q_ph;
                    end
                end
                if ~isempty(src.mapping)
                    qC = src.mapping*qC;
                end
                eq(cells) = eq(cells) - qC;
                src.components{end+1} = qC;

            else
                [eq, src] = addComponentContributions@ReservoirModel(model, cname, eq, component, src, force);
            end
        end

        
        function [state, report] = updateAfterConvergence(model, state0, state, dt, drivingForces)
            [state, report] = updateAfterConvergence@ReservoirModel(model, state0, state, dt, drivingForces);
            if isfield(state, 'wellSol') && isfield(state.wellSol, 'components')
                names = model.getComponentNames();
                for i = 1:numel(names)
                    names{i} = names{i}(isstrprop(names{i}, 'alphanum'));
                end
                ncomp = numel(names);
                for i = 1:numel(state.wellSol)
                    for j = 1:ncomp
                        if state.wellSol(i).status
                            state.wellSol(i).(names{j}) = sum(state.wellSol(i).components(:, j));
                        end
                    end
                end
            end
            
            if isfield(state, 'dz')
                state = rmfield(state, 'dz');
            end
            
            if isfield(state, 'eos')
                state = rmfield(state, 'eos');
            end
        end
        
        function state = setFlag(model, varargin)
            state = model.EOSModel.setFlag(varargin{:});
        end
        
        function [isLiquid, isVapor, is2ph] = getFlag(model, state)
            [isLiquid, isVapor, is2ph] = model.EOSModel.getFlag(state);
        end
        
        function m = PropertyModel(model)
            m = model.EOSModel.PropertyModel;
        end
        
        function scaling = getScalingFactorsCPR(model, problem, names, solver) %#ok
            % Get scaling factors for CPR reduction in `CPRSolverAD`
            %
            % PARAMETERS:
            %   model   - Class instance
            %   problem - `LinearizedProblemAD` which is intended for CPR
            %             preconditioning.
            %   names   - The names of the equations for which the factors are
            %             to be obtained.
            %   solver  - The `LinearSolverAD` class requesting the scaling
            %             factors.
            %
            % RETURNS:
            %   scaling - Cell array with either a scalar scaling factor for
            %             each equation, or a vectogetComponentScalingr of equal length to that 
            %             equation.
            %
            % SEE ALSO
            %   `CPRSolverAD`
            scaling = getScalingFactorsCPR@ReservoirModel(model, problem, names, solver);
%             state = problem.state;
%             rhos = state.rho.*state.s;
%             
%             for i = 1:numel(names)
%                 if strcmpi(names{i}, 'water')
%                     scaling{i} = 1./state.rho(:, 1);
%                 elseif any(strcmpi(names{i}, model.EOSModel.fluid.names))
%                     scaling{i} = 1./max(sum(rhos(:, (1+model.water):end), 2), 1);
%                 end
%             end
        end
    end
    
    methods (Access=protected)
       function [v_comp, tol_comp, names_comp, isComponent] = getComponentConvergenceValues(model, problem)
            % Check convergence criterion of components
            isComponent = false(size(problem.equations));
            cnames = model.getComponentNames();
            names = problem.equationNames;
            
            ncomp = numel(cnames);
            for i = 1:ncomp
                f = strcmpi(names, cnames{i});
                if any(f)
                    isComponent(f) = true;
                    if model.useIncTolComposition
                        names{f} = ['d', cnames{i}];
                    end
                end
            end

            if model.useIncTolComposition
                if problem.iterationNo == 1
                    v_comp = inf(1, ncomp);
                else
                    v_comp = problem.state.dz;
                end
                tol_comp = model.incTolComposition;
            else
                tol_comp = model.nonlinearTolerance;
                v_comp = cellfun(@(x) norm(double(x), inf), problem.equations(isComponent));
            end
            tol_comp = repmat(tol_comp, size(v_comp));
            names_comp = names(isComponent);
       end
       
       function [v_f, tol_f, names_f, isFugacity] = getFugacityConvergenceValues(model, problem)
            % Check fugacity constraints if present
            isFugacity = strcmpi(problem.types, 'fugacity');
            if model.water
                state = problem.state;
                scale = sum(state.s(state.flag == 0, model.water+1:end), 2);
            else
                scale = 1;
            end
            v_f = cellfun(@(x) norm(double(x).*scale, inf), problem.equations(isFugacity));
            tol_f = repmat(model.fugacityTolerance, size(v_f));
            names_f = problem.equationNames(isFugacity);
       end
    end
end


%{
Copyright 2009-2018 SINTEF ICT, Applied Mathematics.

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
