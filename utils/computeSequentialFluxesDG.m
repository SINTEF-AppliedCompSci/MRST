function [q_c, q_f, r_fg, qP_c, qP_f] = computeSequentialFluxesDG(disc, model, state, T, T_all, g, mob, b, sdof, rdof, cdof)
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

    if nargin < 10
        rdof = [];
    end
    if nargin < 11
        cdof = [];
    end
    G = model.G;
    [~, x_f, ~, f] = disc.getCubature((1:G.cells.num)', 'surface');
    flux = sum(state.flux,2);
    % Upstream cells
    [~, ~, e_fv, e_fg] = disc.getSaturationUpwind(f, x_f, T, flux, state, g, mob, sdof, rdof, cdof);
    [~, x_c, e]    = disc.getCubature((1:G.cells.num)', 'volume');
    
    nPh = numel(mob);
    [s_c, sT_c, s_fv, sT_fv, s_fg, sT_fg] = deal(cell(nPh, 1));
    for phNo = 1:nPh
        [s_c{phNo} , sT_c{phNo} ] = disc.evaluateDGVariable(x_c, e           , state, sdof{phNo}, sdof{end});
        [s_fv{phNo}, sT_fv{phNo}] = disc.evaluateDGVariable(x_f, e_fv(:,phNo), state, sdof{phNo}, sdof{end});
        [s_fg{phNo}, sT_fg{phNo}] = disc.evaluateDGVariable(x_f, e_fg(:,phNo), state, sdof{phNo}, sdof{end});
    end
    
    [r_c, r_fv, r_fg] = deal(cell(nPh, 1));
    [r_fv{:}, r_fg{:}] = deal(0);
    if ~isempty(rdof)
        for phNo = 1:nPh
            r_c{phNo}  = disc.evaluateDGVariable(x_c, e           , state, rdof{phNo});
            r_fv{phNo} = disc.evaluateDGVariable(x_f, e_fv(:,phNo), state, rdof{phNo});
            r_fg{phNo} = disc.evaluateDGVariable(x_f, e_fg(:,phNo), state, rdof{phNo});
        end
    end
    
    [r_fv{:}, r_fg{:}] = deal(0);
    if ~isempty(rdof)
        for phNo = 1:nPh
            r_c{phNo}  = disc.evaluateDGVariable(x_c, e           , state, rdof{phNo});
            r_fv{phNo} = disc.evaluateDGVariable(x_f, e_fv(:,phNo), state, rdof{phNo});
            r_fg{phNo} = disc.evaluateDGVariable(x_f, e_fg(:,phNo), state, rdof{phNo});
        end
    end
    [c_c, c_fv, c_fg] = deal(0);
    if ~isempty(cdof)
        c_c  = disc.evaluateDGVariable(x_c, e           , state, cdof);
        c_fv = disc.evaluateDGVariable(x_f, e_fv(:,phNo), state, cdof);
        c_fg = disc.evaluateDGVariable(x_f, e_fg(:,phNo), state, cdof);
    end
    
    [mob_c, mob_fv, mob_fg, b_c, b_fv, b_fg] = deal(cell(nPh, 1));
    sn = @(s, sT) cellfun(@(s, sT) s./sT, s, sT, 'unif', false);
    for phNo = 1:nPh
        mob_c{phNo}  = mob{phNo}(e            , sn(s_c, sT_c)  , r_c{phNo} , c_c);
        mob_fv{phNo} = mob{phNo}(e_fv(:, phNo), sn(s_fv, sT_fv), r_fv{phNo}, c_fv);
        mob_fg{phNo} = mob{phNo}(e_fg(:, phNo), sn(s_fg, sT_fg), r_fg{phNo}, c_fg);
        b_c{phNo}    = b{phNo}(e, s_c{phNo}, r_c{phNo});
        b_fv{phNo}   = b{phNo}(e_fv(:,phNo), s_fv{phNo}, r_fv{phNo});
        b_fg{phNo}   = b{phNo}(e_fg(:,phNo), s_fg{phNo}, r_fg{phNo});
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
        q_cv = b_c{alpha}.*f_c{alpha}.*vT_c(e,:);
        q_fv = b_fv{alpha}.*f_fv{alpha}.*vT(f);
        q_cg = 0;
        q_fg = 0;
        for beta = 1:nPh
            if beta ~= alpha
                q_cg = q_cg + b_c{alpha}.*f_c{alpha}.*mob_c{beta}.*(Kg_c{alpha}(e,:) - Kg_c{beta}(e,:));
                q_fg = q_fg + b_fg{alpha}.*f_fg{alpha}.*mob_fg{beta}.*(Kg_f{alpha}(f) - Kg_f{beta}(f));
            end
        end
        q_c{alpha} = q_cv + q_cg;
        q_f{alpha} = q_fv + q_fg;
    end
    
    [qP_c, qP_f] = deal(0);
    if isprop(model, 'polymer') && model.polymer
        fluid   = model.fluid;
        wIx = model.getPhaseIndex('W');
        a   = fluid.muWMult(fluid.cmax).^(1-fluid.mixPar);
        
        cbar_c  = c_c./fluid.cmax;
        cbar_fv = c_fv./fluid.cmax;
        cbar_fg = c_fg./fluid.cmax;
        
        mobP_c  = mob_c{wIx}./(a + (1-a).*cbar_c);
        mobP_fv = mob_fv{wIx}./(a + (1-a).*cbar_fv);
        mobP_fg = mob_fg{wIx}./(a + (1-a).*cbar_fg);
        
        fP_c  = mobP_c./mobT_c;
        fP_fv = mobP_fv./mobT_fv;
        fP_fg = mobP_fg./mobT_fg;

        qP_cv = b_c{wIx}.*fP_c.*vT_c(e,:);
        qP_fv = b_fv{wIx}.*fP_fv.*vT(f);
        qP_cg = 0;
        qP_fg = 0;
        for beta = 1:nPh
            if beta ~= wIx
                qP_cg = qP_cg + b_c{wIx}.*fP_c.*mob_c{beta}.*(Kg_c{wIx}(e,:) - Kg_c{beta}(e,:));
                qP_fg = qP_fg + b_fg{wIx}.*fP_fg.*mob_fg{beta}.*(Kg_f{wIx}(f) - Kg_f{beta}(f));
            end
        end
        qP_c = (qP_cv + qP_cg).*c_c;
        qP_f = qP_fv.*c_fv + qP_fg.*c_fg;
    end
    
end

% Effective adsorption, depending of desorption or not
function y = effads(c, cmax, model)
   if model.fluid.adsInx == 2
      y = model.fluid.ads(max(c, cmax));
   else
      y = model.fluid.ads(c);
   end
end
%--------------------------------------------------------------------------

% Multipliers due to polymer-----------------------------------------------
function [muWMult, a] = getMobilityMultipliers(model, cMax)
    fluid = model.fluid;
    ads = @(e, c) effads(c, cMax(e), model);
    mixpar = fluid.mixPar;
    cbar   = @(c) c/fluid.cmax;
    a = fluid.muWMult(fluid.cmax).^(1-mixpar);
    b = @(c) 1./(1-cbar(c)+cbar(c)./a);
    % The viscosity multiplier only result from the polymer mixing.
    muWeffMult = @(c) b(c).*fluid.muWMult(c).^mixpar;
    permRed = @(e, c) 1 + ((fluid.rrf-1)./fluid.adsMax).*ads(e, c);
    muWMult = @(e, c) muWeffMult(c).*permRed(e, c);
end
%--------------------------------------------------------------------------
