classdef FacilityModel < PhysicalModel
    properties
        WellModels
        ReservoirModel
    end
    
    methods
        function model = FacilityModel()
            model = model@PhysicalModel([]);
            model.WellModels = {};
        end
        
        function model = setupWells(model, W)
            nw = numel(W);
            if model.getNumberOfWells == 0
                model.WellModels = cell(nw, 1);
                for i = 1:nw
                    % Set up models. SimpleWell for the time being
                    model.WellModels{i} = SimpleWell(W(i));
                end
            else
                assert(model.getNumberOfWells == nw)
                for i = 1:nw
                    % Update with new wells. Typically just a quick
                    % overwrite of existing wells
                    model.WellModels{i}.updateWell(W(i));
                end
            end
        end
        
        function nwell = getNumberOfWells(model)
            nwell = numel(model.WellModels);
        end
        
        function [variables, names, wellmap] = getFacilityPrimaryVariables(model, wellSol)
            nw = model.getNumberOfWells();
            for i = 1:nw
                [v, n] = model.WellModels{i}.getWellPrimaryVariables(wellSol, model.ReservoirModel);
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
            model.ReservoirModel = resModel;
        end
        
        function [srcMass, srcVol, eqs, enames, etypes, wellSol] = getWellContributions(model, wellSol, wellvars, wellMap, p, sat, rho, comp, iteration)
            if isnan(iteration) || iteration < 0
                warning(['Iteration number is not passed on to well model,', ...
                         'this may indicate wellbore pressure-drop will never be updated']);
            end
            
            nw = model.getNumberOfWells();
            for i = 1:nw
                wm = model.WellModels{i};
                
                W = wm.W;
                wc = W.cells;
                pw = p(wc);
                satw = getCellSubset(sat, wc);
                rhow = getCellSubset(rho, wc);
                compw = getComponentCellSubset(comp, wc);
                varw = getVariableSubsetWell(wellvars, wellMap, i);
                % Update pressure
                if iteration == 1
                    wellSol(i) = wm.updateConnectionPressureDrop(wellSol(i), model.ReservoirModel, varw, pw, satw, rhow, compw);
                end
                % Update limits
                wellSol(i) = wm.updateConnectionPressureDrop(wellSol(i), model.ReservoirModel, varw, pw, satw, rhow, compw);
                
                % Set up equations
            end
        end
        
        function wellSol = updateWellSol(model, wellSol)
            
        end
        
        function wc = getWellCells(model)
            
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