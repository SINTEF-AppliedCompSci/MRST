function [s_min, s_max] = getMinMaxPhaseSaturationsFromRelPerm(model, tol, cellInx)
    if nargin < 2
        tol = 1e-6;
    end
    

    nph = sum(model.getActivePhases());
    s = (0:tol:1)';
    
    if nargin < 3
        cellInx = ones(size(s));
    else
        cellInx = repmat(cellInx, size(s));
    end
    
    s_min = zeros(1, nph);
    s_max = ones(1, nph);
    
    if model.water
        krw = model.fluid.krW(s, 'cellInx', cellInx);
        
        ix = model.getPhaseIndex('W');
        s_min(ix) = getMinSat(s, krw);
        s_max(ix) = getMaxSat(s, krw);
    end
    
    if model.oil
        if isfield(model.fluid, 'krO')
            kro = model.fluid.krO(s, 'cellInx', cellInx);
        else
            if model.water
                krow = model.fluid.krOW(s, 'cellInx', cellInx);
            else
                krow = 1;
            end
            
            if model.gas
                if isfield(model.fluid, 'krO')
                    krog = model.fluid.krO(s, 'cellInx', cellInx);
                else
                    krog = model.fluid.krOG(s, 'cellInx', cellInx);
                end
            else
                krog = 1;
            end
            kro = min(krow, krog);
        end
        
        ix = model.getPhaseIndex('O');
        s_min(ix) = getMinSat(s, kro);
        s_max(ix) = getMaxSat(s, kro);
    end
    
    if model.gas
        krg = model.fluid.krG(s, 'cellInx', cellInx);
        
        ix = model.getPhaseIndex('G');        
        s_min(ix) = getMinSat(s, krg);
        s_max(ix) = getMaxSat(s, krg);
    end
end

function s_min = getMinSat(s, kr)
    sub = find(kr == 0, 1, 'last');
    s_min = s(sub);
end

function s_max = getMaxSat(s, kr)
    sub = find(kr == kr(end), 1, 'first');
    s_max = s(sub);
end