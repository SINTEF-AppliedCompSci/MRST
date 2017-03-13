function packed = packPerforationProperties(W, p, mob, rho, dissolved, comp, wellvars, wellvarnames, varmaps, wellmap, ix)
    wc = W.cells;%                          W, p, mob, rho, dissolved, comp, wellvars, addedVars, varmaps, wellMap, i
    packed = struct();
    packed.pressure = p(wc);
    packed.mob = getCellSubset(mob, wc);
    packed.rho = getCellSubset(rho, wc);
    packed.dissolved = getComponentCellSubset(dissolved, wc);
    packed.components = getCellSubset(comp, wc);
    
    % Extra variables outside of standard subset
    varw = getVariableSubsetWell(wellvars, varmaps, ix);
    renum = wellmap(ix, wellmap(ix, :) > 0);
    varw = varw(renum);
    packed.extravars = varw;
    packed.extravars_names = wellvarnames;
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