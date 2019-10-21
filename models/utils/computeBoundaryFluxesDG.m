function [q, rhogKdz] = computeBoundaryFluxesDG(model, state, bc)

    disc = model.discretization;
    G    = model.G;
    % Compute total flux
    flux = sum(state.flux,2);
    vT   = flux./G.faces.areas;
    % Boundary faces
    faces = bc.face;
    % Saturation outside boundary
    s_l   = bc.sat;
    % Face cubature foordinates
    [~, x_bc, ~, face_bc] = disc.getCubature(faces, 'face');
    c_bc = sum(G.faces.neighbors(face_bc,:),2);
    % Mapping from BC faces to global faces
    globFace2BCface        = nan(G.faces.num,1);
    globFace2BCface(faces) = 1:numel(faces);        
    locFaceNo = globFace2BCface(face_bc);
    % Determine injection boundaries
    sgn   = 1 - 2*(G.faces.neighbors(faces, 1) == 0);
    isInj = vT(faces) > 0 & sgn < 0;
    nPh = nnz(model.getActivePhases());
    
    hasOutsideMob = isfield(bc, 'mob');
    hasOutsideRho = isfield(bc, 'rho');
    
    [state_r, state_l] = deal(state);
    
    state_l.s = bc.sat(locFaceNo,:);
    isP = reshape(strcmpi(bc.type, 'pressure'), [], 1);
    if any(isP)
        error('');
        val = bc.value(locFaceNo);
        state_l.pressure = val(isP);
    end
    
    mob_r = model.getProp(state_r, 'Mobility');
    if hasOutsideMob
        error('');
        mob_l = bc.mob(locFaceNo,:);
    else
        mob_l = model.getProp(state_l, 'Mobility');
    end
    
    mob = cell(1,nPh);
    mobT = 0;
    for i = 1:nPh
        mob{i} = mob_l{i}.*isInj + mob_r{i}.*(~isInj);
        mobT   = mobT + mob{i};
    end
    
    rho_r = model.getProp(state_r, 'Density');
    if hasOutsideRho
        error('');
        rho_l = bc.mob(locFaceNo,:);
    else
        rho_l = model.getProp(state_l, 'Density');
    end
    
    rho = cell(1, nPh);
    for i = 1:nPh
        rho{i} = (rho_l{i} + rho_r{i})/2;
    end
    
    
    
%     [s_bc, mob_bc] = deal(cell(nPh, 1));
%     sT_bc = 0;
%     % Evaluate saturations
%     for phNo = 1:nPh
%         % Saturation on inside of boundary
%         s_r = disc.evaluateDGVariable(x_bc, c_bc, state, sdof{phNo});
%         % Upstream saturation
%         s_bc{phNo}   = s_l(locFaceNo,phNo).*isInj(locFaceNo) + s_r.*(~isInj(locFaceNo));
%         % Total saturation
%         sT_bc = sT_bc + s_bc{phNo};
%     end
%     s_bc = cellfun(@(s) s./sT_bc, s_bc, 'unif', false);
    % Compute mobilities
%     mobT_bc = 0;
%     for phNo = 1:nPh
%         mob_bc{phNo} = mob{phNo}(c_bc, s_bc);
%         mobT_bc = mobT_bc + mob_bc{phNo};
%     end
    % Fractional flow functions

    % Update gravity fluxes
    dz = G.cells.centroids(sum(G.faces.neighbors(bc.face,:),2), :) - G.faces.centroids(bc.face,:);
    gvec = model.getGravityVector();
    b_bc = cell(nPh,1);
%     c    = sum(G.faces.neighbors(bc.face,:),2);
%     x    = G.faces.centroids(bc.face,:);
    rhogKdz = cell(1, nPh);
    T = model.operators.T_all(face_bc);
    for i = 1:nPh
        rhogKdz{i} = T.*rho{i}.*(dz*gvec');
    end
    
%     [Kg_bc] = deal(cell(nPh,1));
%     for phNo = 1:nPh
%         Kg_bc{phNo} = T_all.*rhogKdz{phNo};
% %         Kg_c{phNo} = disc.velocityInterp.faceFlux2cellVelocity(Kg_f{phNo});
%         Kg_bc{phNo} = Kg_bc{phNo}./G.faces.areas;
%     end
    q = cell(1,nPh);
    f = cellfun(@(mob) mob./mobT, mob, 'unif', false);
    for alpha = 1:nPh
%         q_vc = b_c{alpha}.*f_bc{alpha}.*vT_c(c,:);
        q_vf = rho{alpha}.*f{alpha}.*vT(face_bc);
%         q_gc = 0;
        q_gf = 0;
        for beta = 1:nPh
            if beta ~= alpha
%                 q_gc = q_gc + b_c{alpha}.*f_c{alpha}.*mob_c{beta}.*(Kg_c{alpha}(c,:) - Kg_c{beta}(c,:));
                q_gf = q_gf + rho{alpha}.*f{alpha}.*mob{beta}.*(rhogKdz{alpha}- rhogKdz{beta});
            end
        end
        q{alpha} = q_vf + q_gf;
    end

end