classdef FacilityModel < PhysicalModel
    properties
        WellModels
        ReservoirModel
    end
    
    properties (Access = protected)
        addedPrimaryVarNames = {};
    end
    methods
        function model = FacilityModel(reservoirModel)
            model = model@PhysicalModel([]);
            model.WellModels = {};
            model.ReservoirModel = reservoirModel;
        end
        
        function model = setupWells(model, W)
            % Set up wells for changed controls or first simulation step
            nw = numel(W);
            if model.getNumberOfWells == 0
                % First time setup
                pvars = cell(nw, 1);
                model.WellModels = cell(nw, 1);
                for i = 1:nw
                    % Set up models. SimpleWell for the time being
                    model.WellModels{i} = SimpleWell(W(i));
                    pvars{i} = model.WellModels{i}.getExtraPrimaryVariableNames();
                end
                model.addedPrimaryVarNames = unique([pvars{:}]);
            else
                assert(model.getNumberOfWells == nw)
                for i = 1:nw
                    % Update with new wells. Typically just a quick
                    % overwrite of existing wells
                    model.WellModels{i}.updateWell(W(i));
                end
                % Update wellSol as well
            end
        end
        
        function nwell = getNumberOfWells(model)
            % Compute number of wells in facility
            nwell = numel(model.WellModels);
        end
        
        function names = getPrimaryVariableNames(model, wellSol)
            % This includes both the basic variables, and the variables
            % added by complex wells (if any)
            names = [model.getBasicPrimaryVariableNames(), model.addedPrimaryVarNames];
        end
        
        function names = getBasicPrimaryVariableNames(model)
            % Basic primary variables are phase rates + bhp for active
            % phases in the model.
            actPh = model.ReservoirModel.getActivePhases();
            names = {'qWs', 'qOs', 'qGs', 'bhp'};
            names = names([actPh, true]);
        end

        function [rates, bhp, names] = getBasicPrimaryVariables(model, wellSol)
            actPh = model.ReservoirModel.getActivePhases();
            bhp = vertcat(wellSol.bhp);
            qWs = vertcat(wellSol.qWs);
            qOs = vertcat(wellSol.qOs);
            qGs = vertcat(wellSol.qGs);
            rates = {qWs, qOs, qGs};
            rates = rates(actPh);
            
            names = model.getBasicPrimaryVariableNames();
        end
        
        function [vars, names, wellmap] = getExtraPrimaryVariables(model, wellSol)
            % Extra primary variables are variables required by more
            % advanced wells that are in addition to the basic facility
            % variables (rates + bhp).
            nw = model.getNumberOfWells();
            for i = 1:nw
                [v, n] = model.WellModels{i}.getExtraPrimaryVariables(wellSol, model.ReservoirModel);
                if i == 1
                    % Currently assuming that all wells are of the same
                    % type
                    nv = numel(v);
                    vars = cell(nw, nv);
                    names = n;
                end
                for j = 1:numel(v)
                    vars{i, j} = v{j};
                end

            end
            
            variables = cell(1, nv);
            wellmap = cell(1, nv);
            for j = 1:nv
                variables{j} = vertcat(vars{:, j});
                wellmap{j} = rldecode((1:nw)', cellfun(@numel, vars(:, j)));
            end
        end
        
        function model = setReservoirModel(model, resModel)
            % Store pointer to reservoir model (used to figure out which
            % components and phases are present)
            model.ReservoirModel = resModel;
        end
        
        function [srcMass, srcVol, eqs, ctrleq, enames, etypes, wellSol] = getWellContributions(model, wellSol, qWell, bhp, wellvars, wellMap, p, mob, rho, comp, iteration)
            % Get the source terms due to the wells, control and well
            % equations and updated well sol. Main gateway for adding wells
            % to a set of equations.
            if isnan(iteration) || iteration < 0
                warning(['Iteration number is not passed on to well model,', ...
                         'this may indicate wellbore pressure-drop will never be updated']);
            end
            
            nw = model.getNumberOfWells();
            
            allEqs = cell(nw, 1);
            allCtrl = cell(nw, 1);
            
            allVol = cell(nw, 1);
            allMass = cell(nw, 1);
            for i = 1:nw
                wm = model.WellModels{i};
                [enames, etypes] = wm.getWellEquationNames(model.ReservoirModel);
                
                W = wm.W;
                wc = W.cells;
                pw = p(wc);
                mobw = getCellSubset(mob, wc);
                rhow = getCellSubset(rho, wc);
                compw = getComponentCellSubset(comp, wc);
                varw = getVariableSubsetWell(wellvars, wellMap, i);
                qw = cellfun(@(x) x(i), qWell, 'uniformoutput', false);
                bh = bhp(i);
                % Update pressure
                if iteration == 1
                    wellSol(i) = wm.updateConnectionPressureDrop(wellSol(i), model.ReservoirModel, qw, bh, varw, pw, mobw, rhow, compw);
                end
                % Update limits
                [qw, bh, wellSol(i), ok] = wm.updateLimits(wellSol(i), model.ReservoirModel, qw, bh, varw, pw, mobw, rhow, compw);
                if ~ok
                    bhp(i) = bh;
                    for phNo = 1:numel(qw)
                        qWell{phNo}(i) = qw{phNo};
                    end
                end
               % Set up well equations and source terms
               [allEqs{i}, allCtrl{i}, allMass{i}, allVol{i}, wellSol(i)] =...
                   wm.computeWellEquations(wellSol(i), model.ReservoirModel, qw, bh, varw, pw, mobw, rhow, compw);
            end
            nPh = nnz(model.ReservoirModel.getActivePhases);
            [srcMass, srcVol, eqs] = deal(cell(1, nPh));
            for phNo = 1:nPh
                srcMass{phNo} = combineCellData(allMass, phNo);
                srcVol{phNo} = combineCellData(allVol, phNo);
                eqs{phNo} = combineCellData(allEqs, phNo);
            end
            ctrleq = vertcat(allCtrl{:});
        end
        

        function wellSol = updateWellSolAfterStep(model, resmodel, wellSol)
            % Figure out if wells are shut, or changed ontrols
            for wno = 1:numel(wellSol)
                wm = model.WellModels{wno};
                wellSol(wno) = wm.updateWellSolAfterStep(resmodel, wellSol(wno));
            end
        end
        
        function wc = getWellCells(model)
            c = cellfun(@(x) x.W.cells, model.WellModels, 'UniformOutput', false);
            wc = vertcat(c{:});
        end
        
    end
end

function celldata = getComponentCellSubset(celldata, wc)
    for i = 1:numel(celldata)
        for j = 1:numel(celldata{i});
            if ~isempty(celldata{i}{j})
                celldata{i}{j} = celldata{i}{j}(wc);
            end
        end
    end
end

function d = combineCellData(data, ix)
    d = cellfun(@(x) x{ix}, data, 'UniformOutput', false);
    d = vertcat(d{:});
end

function subset = getCellSubset(celldata, wc)
    subset = cell(size(celldata));
    for i = 1:numel(subset)
        if ~isempty(celldata{i})
            subset{i} = celldata{i}(wc);
        end
    end
end

function subset = getVariableSubsetWell(vars, wellMap, ix)
    subset = cell(size(vars));
    for i = 1:numel(subset)
        subset{i} = vars{i}(wellMap{i} == ix);
    end
end