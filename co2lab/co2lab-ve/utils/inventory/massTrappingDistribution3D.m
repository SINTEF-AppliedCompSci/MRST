function masses = massTrappingDistribution3D(s, smax, pressure, rs, G, ...
                                             fluid, rock, res_water, res_gas, varargin)
% Compute the trapping distribution of CO2 in each cell of a 3D gird
%
% SYNOPSIS:
%   function masses = massTrappingDistribution3D(s, smax, pressure, rs, Gt, ...
%                                                fluidVE, rockVE, res_water, ...
%                                                res_gas, varargin)
%
% DESCRIPTION:
%
% PARAMETERS:
%   s          - current CO2 saturation per cell
%   smax       - historically maximum OC2 saturation per cell
%   pressure   - cell-wise pressure
%   rs         - fraction of dissolved CO2 in brine per cell
%   G          - 3D grid
%   fluid      - fluid object
%   rock       - rock object
%   res_water  - max value for residual water saturation
%   res_gas    - max value for residual CO2 saturation 
%   optional values:
%      dh         - subtrapping capacity (support not yet implemented)
%      trapstruct - trapping structure (if not provided, it will be calculated)
% 
% RETURNS:
%   masses - vector with 7 components, representing:
%            masses[1] : mass of dissolved gas, per cell
%            masses[2] : mass of gas that is both structurally and residually trapped
%            masses[3] : mass of gas that is residually (but not structurally) trapped
%            masses[4] : mass of non-trapped gas that will be residually trapped
%            masses[5] : mass of structurally trapped gas, not counting the gas that 
%                        will eventually be residually trapped
%            masses[6] : mass of subscale trapped gas (if 'dh' is nonempty)
%            masses[7] : mass of 'free' gas (i.e. not trapped in any way)

%{
Copyright 2009-2025 SINTEF Digital, Mathematics & Cybernetics.

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

    
    opt.dh = 0;
    opt.trapstruct = [];
    opt = merge_options(opt, varargin{:});

    % warning if dh ~= 0, since support has not yet been implemented
    if opt.dh ~= 0
        warning('dh option not yet supported in massTrappingDistribution3DNew');
    end
    
    Gt = topSurfaceGrid(G);
    
    krG = get_imb_relperm_curve(fluid, 'krG');
    krW = get_imb_relperm_curve(fluid, 'krW');
    
    % srg_max = find_rgmax(krG);
    % srw_max = find_rgmax(krW);
    
    pvMult = 1;
    if isfield(fluid, 'pvMultR')
        pvMult = fluid.pvMultR(pressure);
    end
    
    pv = rock.poro .* G.cells.volumes .* pvMult;
    rhoG = fluid.rhoGS .* fluid.bG(pressure);
    mfac = pv .* rhoG;
    
    % % compute interface depths for each pillar
    h_max = find_lowest_cell_z(G, Gt, smax > 1e-2);
    h     = find_lowest_cell_z(G, Gt, krG(s) > 1e-3);
    
    % find trapcells
    trapcells = identify_3D_trapcells(Gt, opt.trapstruct);
    
    % compute mobile and immobile saturation
    s_immob = min(land_immobile_saturation(smax, res_gas, res_water), s);
    s_mob = s - s_immob;
    
    total_undissolved_mass = sum(mfac .* s); % used for sanity check
    
    %% --- computing masses to go in reporting ---
    
    % dissolved CO2
    disgas_mass = sum(fluid.rhoGS .* fluid.bW(pressure) .* pv .* rs .* (1-s));
    
    % immobilized, structurally trapped CO2
    immob_struct_mass = sum(mfac .* s_immob .* trapcells);
    
    % residual CO2 (below flowing plume, i.e. between h and hmax
    h_cells = pillar_values_to_cells(Gt, h);
    
    topz_cells = pillar_values_to_cells(Gt, Gt.cells.z);
    d = G.cells.centroids(:,3) - topz_cells;
    
    res_cells = d > h_cells; % cells that are below the moving plume and not in a trap
    residual_mass = sum(mfac .* s_immob .* (res_cells & ~trapcells)); 
    
    % immobilized, free CO2 (other than residual, i.e. inside moving plume)
    immob_free_mass = sum(mfac .* s_immob .* (~res_cells & ~trapcells));
    
    % mobile structural trapped CO2
    mobile_struct_mass = sum(mfac .* s_mob .* trapcells);
    
    % CO2 in subtraps
    subtrap_mass = 0; % @@ not implemented yet
    
    % mobile, free-flowing CO2
    mobile_free_mass = sum(mfac .* s_mob .* ~trapcells);

    % assemble return variable
    masses = max([value(disgas_mass), ...         
                  value(immob_struct_mass), ...   
                  value(residual_mass), ...
                  value(immob_free_mass), ... 
                  value(mobile_struct_mass), ...  
                  value(subtrap_mass), ...        
                  value(mobile_free_mass)], 0);   

    % sanity check
    if(abs(sum(masses(2:end))-total_undissolved_mass) > 1e-3 * total_undissolved_mass)
        abs(sum(masses(2:end))-total_undissolved_mass)
        abs(sum(masses(2:end))-total_undissolved_mass) - 1e-3 * total_undissolved_mass
        disp('There is a mismatch between mass calculations');
    end
    
end

function cvals = pillar_values_to_cells(Gt, pvals)

    cmap = rldecode((1:Gt.cells.num)', diff(Gt.cells.columnPos));
    cvals(Gt.columns.cells) = pvals(cmap);
    cvals = cvals(:);
end


function tc_ind = identify_3D_trapcells(Gt, trapstruct)
    
    if isempty(trapstruct)
        trapstruct = trapAnalysis(Gt, false);
    end
    
    spill = Gt.cells.z;
    ix = find(trapstruct.traps > 0);
    spill(ix) = trapstruct.trap_z(trapstruct.traps(ix));
    
    field = pillar_values_to_cells(Gt, spill);
    tc_ind = field(:) >= Gt.parent.cells.centroids(:,3);
    
end

    
function krfun = get_imb_relperm_curve(fluid, field)
    
    krfun = fluid.(field);
    if iscell(krfun)
        assert(isfield(fluid, 'krHyst')); % should only happen for fluid models
                                          % with hysteresis
        krfun = krfun{end};
    end
end

% function rg_max = find_rgmax(krfun)
% % @@ at the moment, a quick, dirty and imprecise implementation, but 
% % probably sufficient?
%     s = linspace(0, 1, 3000);
%     kr = krfun(s);
%     ix = find(kr < sqrt(eps), 1, 'last');
%     rg_max = s(ix);
% end


function s_immob = land_immobile_saturation(smax, srg_max, srw_max)

    C = 1 / (srg_max+eps) - 1 ./ (1 - srw_max);
    
    s_immob = smax ./ (1 + C * smax);
    
end



function h_lowest = find_lowest_cell_z(G, Gt, ind)
    
    cmap = rldecode((1:Gt.cells.num)', diff(Gt.cells.columnPos));
    
    colind(Gt.columns.cells) = cmap;
    zvals = G.cells.centroids(:,3);
    
    mat = [colind(:), zvals]; % first column is index of vertical pillar, second the
                              % cell z-values
    
    mat(~ind, :) = []; % remove cells that do not fulfill given criterion
    mat = sortrows(mat, [1, 2]);
    remove = find(mat(1:end-1, 1) == mat(2:end, 1));
    mat(remove, :) = []; % keep only last entry (with largest z value)
    
    h_lowest = Gt.cells.z;
    h_lowest(mat(:, 1)) = mat(:, 2); 
    h_lowest = h_lowest - Gt.cells.z;
end

