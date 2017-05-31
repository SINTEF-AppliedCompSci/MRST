function [src, bc] = computeSourcesAndBoundaryConditionsAD(model, pressure, s, mob, rho, dissolved, forces)

    if isempty(forces.bc) && isempty(forces.src)
        return
    end
    rhoS = model.getSurfaceDensities();
    b = phaseDensitiesTobfactor(rho, rhoS, dissolved);
    
    [qVolBC, qVolSRC] = deal(cell(1, numel(rhoS)));
    [bcCells, srcCells] = deal([]);
    
    if ~isempty(forces.bc)
        % Setup the fluxes from the boundary condition
        [qVolBC, BCTocellMap, bcCells] = getBoundaryConditionFluxesAD(model, pressure, b, mob, s, forces.bc);
                
        qVolBC = fixRepeats(qVolBC, BCTocellMap, bcCells);
    end
    
    if ~isempty(forces.src)
        % Fluxes from source terms
        [qVolSRC, srcCells] = getSourceFluxesAD(model, mob, s, forces.src);
    end
    src = getContributionsStruct(model, forces.src, qVolSRC, b, rhoS, srcCells, dissolved);
    bc = getContributionsStruct(model, forces.bc, qVolBC, b, rhoS, bcCells, dissolved);
end

function val = fixRepeats(val, map, cells)
    for i = 1:numel(val)
        v = map*val{i};
        val{i} = v(cells);
    end
end

function src = getContributionsStruct(model, force, q_s, b, rhoS, cells, dissolved)
    nPh = numel(q_s);
    q_r = q_s;
    for i = 1:nPh
        q_r{i} = q_r{i}./b{i}(cells);
    end
    
    if ~isempty(dissolved)
        assert(0);
    end
    
    srcMass = q_s;
    for i = 1:nPh
        srcMass{i} = srcMass{i}.*rhoS(i);
    end
    src = struct('phaseMass',   {srcMass}, ...
                     'phaseVolume', {q_r}, ...
                     'components',  {[]}, ...
                     'sourceCells', cells);
end