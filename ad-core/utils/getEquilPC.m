function [PC, pc_sign, pc_scalers] = getEquilPC(model, satnum, cells)

    actPh = model.getActivePhases();
    nPh = sum(actPh);
    f = model.fluid;

    nc = numel(cells);
    PC = cell(1, nPh);
    pc_sign = ones(1, nPh);
    pc_scalers = ones(nc, nPh);

    pcimpl = model.FlowPropertyFunctions.CapillaryPressure;
    
    if model.water
        ix = model.getPhaseIndex('W');
        
        pc_sign(ix) = -1;
        if isfield(model.fluid, 'pcOW')
            pcOW = getFunction(f, 'pcOW', satnum);
            PC{ix} = @(S) pcOW(S);
            if pcimpl.hasJFunctionScaler('OW')
                ratio = pcimpl.getJFunctionStaticRatio(model);
                pc_scalers(:, ix) = ratio(cells).*pcimpl.getSurfaceTension('OW');
            end
        else
            PC{ix} = @(S) 0*S;
        end
    end
    if model.oil
        ix = model.getPhaseIndex('O');
        PC{ix} = @(S) 0*S;
    end
    if model.gas
        ix = model.getPhaseIndex('G');
        pc_sign(ix) = 1;
        if isfield(model.fluid, 'pcOG')
            pcOG = getFunction(f, 'pcOG', satnum);
            PC{ix} = @(S) pcOG(S);
            if pcimpl.hasJFunctionScaler('OG')
                ratio = pcimpl.getJFunctionStaticRatio(model);
                pc_scalers(:, ix) = ratio(cells).*pcimpl.getSurfaceTension('OG');
            end
        elseif ~model.oil && isfield(model.fluid, 'pcWG')
            pcWG = getFunction(f, 'pcWG', satnum);
            PC{ix} = @(S) pcWG(S);
        else
            PC{ix} = @(S) 0*S;
        end
    end
end

function f = getFunction(fluid, fld, reg)
    f = fluid.(fld);
    if iscell(f)
        f = f{reg};
    end 
end
