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
    
    if ~isempty(dissolved) && ~isempty(force)
        q_s0 = q_s;
        q_t = 0;
        for i = 1:numel(q_s)
            q_t = q_t + double(q_s{i});
        end
        isInj = q_t > 0;
        for i = 1:nPh
            for j = 1:nPh
                % Add dissolution of component i into phase j
                r_cell = dissolved{i}{j};
                if isempty(r_cell)
                    continue
                end
                
                if isfield(force, 'dissolution')
                    ds_mat = force.dissolution;
                    if ndims(ds_mat) == 3
                        % One value per cell, per phase and component
                        r_inj = ds_mat(:, i, j);
                    else
                        % One value for all cells, per phase and component
                        r_inj = ds_mat(i, j);
                    end
                else
                    r_inj = 0;
                end
                r = isInj.*r_inj + ~isInj.*r_cell(cells);
                q_s{i} = q_s{i} + r.*q_s0{j};
            end
        end
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