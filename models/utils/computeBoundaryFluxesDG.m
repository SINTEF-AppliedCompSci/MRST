function [q_bc, g] = computeBoundaryFluxesDG(model, state, bc, T_all, g, mob, b, rho, sdof, rdof)

    disc = model.disc;
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
    nPh = numel(sdof)-1;
    [s_bc, mob_bc] = deal(cell(nPh, 1));
    sT_bc = 0;
    % Evaluate saturations
    for phNo = 1:nPh
        % Saturation on inside of boundary
        s_r = disc.evaluateDGVariable(x_bc, c_bc, state, sdof{phNo});
        % Upstream saturation
        s_bc{phNo}   = s_l(locFaceNo,phNo).*isInj(locFaceNo) + s_r.*(~isInj(locFaceNo));
        % Total saturation
        sT_bc = sT_bc + s_bc{phNo};
    end
    s_bc = cellfun(@(s) s./sT_bc, s_bc, 'unif', false);
    % Compute mobilities
    mobT_bc = 0;
    for phNo = 1:nPh
        mob_bc{phNo} = mob{phNo}(c_bc, s_bc);
        mobT_bc = mobT_bc + mob_bc{phNo};
    end
    % Fractional flow functions
    f_bc = cellfun(@(mob) mob./mobT_bc, mob_bc, 'unif', false);
    % Update gravity fluxes
    dz = G.cells.centroids(sum(G.faces.neighbors(bc.face,:),2), :) - G.faces.centroids(bc.face,:);
    gvec = model.getGravityVector();
    b_bc = cell(nPh,1);
    c    = sum(G.faces.neighbors(bc.face,:),2);
    x    = G.faces.centroids(bc.face,:);
    for phNo = 1:nPh
        b_bc{phNo} = b{phNo}(c_bc, s_bc{phNo});
        
        s      = disc.evaluateDGVariable(x, c, state, sdof{phNo});
        rho_bc = rho{phNo}(c, s);
        g{phNo}(bc.face) = rho_bc.*(dz*gvec');
    end
    
    [Kg_bc] = deal(cell(nPh,1));
    for phNo = 1:nPh
        Kg_bc{phNo} = T_all.*g{phNo};
%         Kg_c{phNo} = disc.velocityInterp.faceFlux2cellVelocity(Kg_f{phNo});
        Kg_bc{phNo} = Kg_bc{phNo}./G.faces.areas;
    end
    q_bc = deal(cell(nPh,1));
    for alpha = 1:nPh
%         q_vc = b_c{alpha}.*f_bc{alpha}.*vT_c(c,:);
        q_vf = b_bc{alpha}.*f_bc{alpha}.*vT(face_bc);
%         q_gc = 0;
        q_gf = 0;
        for beta = 1:nPh
            if beta ~= alpha
%                 q_gc = q_gc + b_c{alpha}.*f_c{alpha}.*mob_c{beta}.*(Kg_c{alpha}(c,:) - Kg_c{beta}(c,:));
                q_gf = q_gf + b_bc{alpha}.*f_bc{alpha}.*mob_bc{beta}.*(Kg_bc{alpha}(face_bc) - Kg_bc{beta}(face_bc));
            end
        end
        q_bc{alpha} = q_vf + q_gf;
    end

end