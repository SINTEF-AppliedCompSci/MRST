function [q, rhogKdz] = computeBoundaryFluxesDG(model, state, bc)
%Undocumented Utility Function

%{
Copyright 2009-2022 SINTEF Digital, Mathematics & Cybernetics.

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

    disc = model.operators.discretization;
    G    = model.G;
    % Compute total flux
    
    vT = sum(state.flux,2);
    
%     vT   = bc.value;
    % Boundary faces
    faces = bc.face;
    % Face cubature foordinates
    [~, ~, ~, faceNo] = disc.getCubature(faces, 'face');
    % Mapping from BC faces to global faces
    globFace2BCface        = nan(G.faces.num,1);
    globFace2BCface(faces) = 1:numel(faces);        
    locFaceNo = globFace2BCface(faceNo);
    % Determine injection boundaries
    sgn   = 2*(G.faces.neighbors(faces, 1) == 0) - 1;
    isInj = vT(faces) > 0 & sgn > 0;
    isInj = isInj(locFaceNo);
    nPh = nnz(model.getActivePhases());
    
    hasOutsideMob = isfield(bc, 'mob');
    hasOutsideRho = isfield(bc, 'rho');
    
    [state_r, state_l] = deal(state);
    
    isP = reshape(strcmpi(bc.type, 'pressure'), [], 1);
    if any(isP)
        error('');
        val = bc.value(locFaceNo);
        state_l.pressure = val(isP);
    end
    
    mob_r = model.getProp(state_r, 'Mobility');
    state_l = model.reduceState(state_l);
    state_l.s = bc.sat(locFaceNo,:);
    if hasOutsideMob
        error('');
        mob_l = bc.mob(locFaceNo,:);
    else
        mob_l = model.getProp(state_l, 'Mobility');
    end
   
    rho_r = model.getProp(state_r, 'Density');
    if hasOutsideRho
        error('');
        rho_l = bc.rho(locFaceNo,:);
    else
        rho_l = model.getProp(state_l, 'Density');
    end
    
    [mob, rho] = deal(cell(1,nPh));
    mobT = 0;
    upw = @(vl, vr) vl.*isInj + vr.*(~isInj);
    for i = 1:nPh
        mob{i} = upw(mob_l{i}, mob_r{i});
        rho{i} = upw(rho_l{i}, rho_r{i});
        mobT   = mobT + mob{i};
    end

    % Update gravity fluxes
    dz = G.cells.centroids(sum(G.faces.neighbors(faceNo,:),2), :) - G.faces.centroids(faceNo,:);
    gvec = model.getGravityVector();
    rhogKdz = cell(1, nPh);
    T = model.operators.T_all(faceNo);
    for i = 1:nPh
        rhogKdz{i} = T.*rho{i}.*(dz*gvec');
    end
    
%     sT = model.getProp(state_l, 'TotalSaturation');
    
    q = cell(1,nPh);
    f = cellfun(@(mob) mob./mobT, mob, 'unif', false);
    for alpha = 1:nPh
        qV = rho{alpha}.*f{alpha}.*vT(faceNo);
        qG = 0;
        for beta = 1:nPh
            if beta ~= alpha
                qG = qG + rho{alpha}.*f{alpha}.*mob{beta}.*(rhogKdz{alpha}- rhogKdz{beta});
            end
        end
        q{alpha} = qV + qG;
    end

end
