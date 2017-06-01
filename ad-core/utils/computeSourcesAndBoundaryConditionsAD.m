function [src, bc] = computeSourcesAndBoundaryConditionsAD(model, pressure, s, mob, rho, dissolved, forces)
    hasBC  = ~isempty(forces.bc);
    hasSRC = ~isempty(forces.src);
    rhoS = model.getSurfaceDensities();
    if hasBC || hasSRC
        b = phaseDensitiesTobfactor(rho, rhoS, dissolved);
        if hasBC
            % Setup the fluxes from the boundary condition
            [qVolBC, BCTocellMap, bcCells] = getBoundaryConditionFluxesAD(model, pressure, b, mob, s, forces.bc);
        end

        if hasSRC
            % Fluxes from source terms
            [qVolSRC, srcCells] = getSourceFluxesAD(model, mob, s, forces.src);
        end
    else
        [qVolBC, qVolSRC, b] = deal(cell(1, numel(mob)));
        [bcCells, srcCells, BCTocellMap] = deal([]);
    end
    src = getContributionsStruct(forces.src, qVolSRC, b, rhoS, srcCells, dissolved);
    bc = getContributionsStruct(forces.bc, qVolBC, b, rhoS, bcCells, dissolved, BCTocellMap);
end

function src = getContributionsStruct(force, q_s, b, rhoS, cells, dissolved, map)
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
                    % Note: If dissolution is specified for rate
                    % sources/bc, the total injected mass will change as
                    % the dissolved mass is not accounted for. This feature
                    % is primarily intended for pressure bc.
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
    if nargin > 6
        mm = map(cells, :);
        for i = 1:numel(q_s)
            q_s{i} = mm*q_s{i};
            q_r{i} = mm*q_r{i};
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