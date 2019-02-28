function [q_c, q_f, r_c, r_fg] = computeSequentialFluxesDG(disc, model, state, T, T_all, g, mob, b, sdof, rdof)

    G = model.G;
    [~, x_f, ~, f] = disc.getCubature((1:G.cells.num)', 'surface');
    flux = sum(state.flux,2);
    % Upstream cells
    [~, ~, c_fv, c_fg] = disc.getSaturationUpwind(f, x_f, T, flux, state, g, mob, sdof, rdof);
    [~, x_c, c]    = disc.getCubature((1:G.cells.num)', 'volume');
    
    nPh = numel(mob);
    [s_c, sT_c, s_fv, sT_fv, s_fg, sT_fg] = deal(cell(nPh, 1));
    for phNo = 1:nPh
        [s_c{phNo} , sT_c{phNo} ] = disc.evaluateDGVariable(x_c, c           , state, sdof{phNo}, sdof{end});
        [s_fv{phNo}, sT_fv{phNo}] = disc.evaluateDGVariable(x_f, c_fv(:,phNo), state, sdof{phNo}, sdof{end});
        [s_fg{phNo}, sT_fg{phNo}] = disc.evaluateDGVariable(x_f, c_fg(:,phNo), state, sdof{phNo}, sdof{end});
    end
    
    [r_c, r_fv, r_fg] = deal(cell(nPh, 1));
    [r_fv{:}, r_fg{:}] = deal(0);
    if ~isempty(rdof)
        for phNo = 1:nPh
            r_c{phNo}  = disc.evaluateDGVariable(x_c, c           , state, rdof{phNo});
            r_fv{phNo} = disc.evaluateDGVariable(x_f, c_fv(:,phNo), state, rdof{phNo});
            r_fg{phNo} = disc.evaluateDGVariable(x_f, c_fg(:,phNo), state, rdof{phNo});
        end
    end
    
    [mob_c, mob_fv, mob_fg, b_c, b_fv, b_fg] = deal(cell(nPh, 1));
    sn = @(s, sT) cellfun(@(s, sT) s./sT, s, sT, 'unif', false);
    for phNo = 1:nPh
        mob_c{phNo}  = mob{phNo}(c            , sn(s_c, sT_c)  , r_c{phNo} );
        mob_fv{phNo} = mob{phNo}(c_fv(:, phNo), sn(s_fv, sT_fv), r_fv{phNo});
        mob_fg{phNo} = mob{phNo}(c_fg(:, phNo), sn(s_fg, sT_fg), r_fg{phNo});
        b_c{phNo}    = b{phNo}(c, s_c{phNo}, r_c{phNo});
        b_fv{phNo}   = b{phNo}(c_fv(:,phNo), s_fv{phNo}, r_fv{phNo});
        b_fg{phNo}   = b{phNo}(c_fg(:,phNo), s_fg{phNo}, r_fg{phNo});
    end
    
    [mobT_c, mobT_fv, mobT_fg] = deal(0);
    for phNo = 1:nPh
        mobT_c  = mobT_c + mob_c{phNo};
        mobT_fv = mobT_fv + mob_fv{phNo};
        mobT_fg = mobT_fg + mob_fg{phNo};
    end
    
    [f_c, f_fv, f_fg] = deal(cell(nPh, 1));
    for phNo = 1:nPh
        f_c{phNo}  = sT_c{phNo}.*mob_c{phNo}./mobT_c;
        f_fv{phNo} = sT_fv{phNo}.*mob_fv{phNo}./mobT_fv;
        f_fg{phNo} = sT_fg{phNo}.*mob_fg{phNo}./mobT_fg;
    end
    
    [Kg_c, Kg_f] = deal(cell(nPh,1));
    for phNo = 1:nPh
        Kg_f{phNo} = T_all.*g{phNo};
        Kg_c{phNo} = disc.velocityInterp.faceFlux2cellVelocity(Kg_f{phNo});
        Kg_f{phNo} = Kg_f{phNo}./G.faces.areas;
    end
    
    vT_c = disc.velocityInterp.faceFlux2cellVelocity(flux);
    vT   = flux./G.faces.areas;
    [q_c, q_f] = deal(cell(nPh,1));
    for alpha = 1:nPh
        q_vc = b_c{alpha}.*f_c{alpha}.*vT_c(c,:);
        q_vf = b_fv{alpha}.*f_fv{alpha}.*vT(f);
        q_gc = 0;
        q_gf = 0;
        for beta = 1:nPh
            if beta ~= alpha
                q_gc = q_gc + b_c{alpha}.*f_c{alpha}.*mob_c{beta}.*(Kg_c{alpha}(c,:) - Kg_c{beta}(c,:));
                q_gf = q_gf + b_fg{alpha}.*f_fg{alpha}.*mob_fg{beta}.*(Kg_f{alpha}(f) - Kg_f{beta}(f));
            end
        end
        q_c{alpha} = q_vc + q_gc;
        q_f{alpha} = q_vf + q_gf;
    end
    
end