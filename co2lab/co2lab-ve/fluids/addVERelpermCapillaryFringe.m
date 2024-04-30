function fluid = addVERelpermCapillaryFringe(fluid, Gt, rock2D, invPc3D, kr3D, varargin)
%Undocumented Utility Function

%{
Copyright 2009-2024 SINTEF Digital, Mathematics & Cybernetics.

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

    % type can be 'linear cap.', 'S table', 'P-scaled table' or 'P-K-scaled table'.
    opt = struct('type', 'P-scaled table', 'samples', 2000);
    opt = merge_options(opt, varargin{:});

    % Warn user if data doesn't fulfill the assumptions for the model he/she
    % has chosen
    issue_warnings(opt.type, fluid, Gt, rock2D);

    % Determine upper bounds for sampled tables
    Pmax = determine_upper_cap_press_limit(fluid, Gt, invPc3D);
    drho_surf = fluid.rhoWS - fluid.rhoGS;
    
    switch opt.type
      case 'linear cap.'
        error('linear cap. unimplemented.'); % @@@@@
      case 'S table'
        table = make_CO2_table_h_based(invPc3D, kr3D, Gt, opt.samples, Pmax, drho_surf);
        fluid.pcWG = @(sg, p, varargin) pcWG_htable(sg, table, fluid, Gt.cells.H, varargin{:});
        fluid.krG = @(sg, p, varargin) krG_htable(sg, table, fluid, Gt.cells.H, varargin{:});
        fluid.krW = @(sw, p, varargin) krW_simple(sw, fluid, varargin{:});
      case 'P-scaled table'
        table = make_CO2_table_p_based(invPc3D, kr3D, opt.samples, Pmax);
        fluid.pcWG = @(sg, p, varargin) pcWG_ptable(sg, p, table, fluid, Gt.cells.H, varargin{:});
        fluid.krG = @(sg, p, varargin) krG_ptable(sg, p, table, fluid, Gt.cells.H, varargin{:});
        fluid.krW = @(sw, p, varargin) krW_simple(sw, fluid, varargin{:});
      case 'P-K-scaled table'
        kscale = sqrt(rock.poro ./ (rock.perm)) * fluid.surface_tension;
        table = make_CO2_table_p_based(invPc3D, kr3D, opt.samples, Pmax / kscale);
        fluid.pcWG = @(sg, p, varargin) ...
            pcWG_ptable(sg, p, table, fluid, Gt.cells.H, 'kscale', kscale, varargin{:});
        fluid.krG = @(sg, p, varargin) ...
            krG_ptable(sg, p, table, fluid, Gt.cells.H, 'kscale', kscale, varargin{:});
        fluid.krW = @(sw, p, varargin) krW_simple(sw, fluid, varargin{:});
      otherwise
        error('Unknown VE relperm model');
    end
    
    % add fine-scale functions to fluid object
    fluid.invPc3D = invPc3D;
    fluid.kr3D = kr3D;
end

% ----------------------------------------------------------------------------
function kr = krW_simple(sw, fluid, varargin)
% For the time being, this is a crude approximation built on a sharp interface
% assumption.  A more accurate approach would use sampled tables for water as well.
        
    opt = merge_options(struct('sGmax', []), varargin{:});
    if isempty(opt.sGmax)
        opt.sGmax = 1-sw; 
    end
        
    sg = free_sg(1-sw, opt.sGmax, fluid.res_water, fluid.res_gas);
    
    sw_eff = sw - (sg./(1-fluid.res_water)) .* fluid.res_water;
    sw_eff(sw_eff<0) = 0; % Should not logically happen, but just in case
    
    kr = sw_eff;
    
end

% ----------------------------------------------------------------------------
function kr = krG_ptable(sg, p, table, fluid, H, varargin)
    opt = merge_options(struct('sGmax', [], 'kscale', 1), varargin{:});
    if isempty(opt.sGmax)
        opt.sGmax = sg;
    end
    
    drho = fluid.rhoW(p) - fluid.rhoG(p);
    sg = free_sg(sg, opt.sGmax, fluid.res_water, fluid.res_gas);
    
    [SP, H_trunc] = compute_untruncated_SP(H, drho, sg, table);
    
    kr = interpTable(table.SP, table.krP, SP ./ opt.kscale);
    
    if any(value(H_trunc))
        % adjust for truncated plumes

        p_trunc = H_trunc .* drho .* norm(gravity);
        kr_trunc = interpTable(table.p, table.krP, p_trunc);
        kr = kr - kr_trunc;
    end
    
    kr = kr ./ (H .* drho .* norm(gravity));
    
end

% ----------------------------------------------------------------------------
function pc = pcWG_ptable(sg, p, table, fluid, H, varargin)
    opt = merge_options(struct('sGmax', [], 'kscale', 1), varargin{:});    
    if isempty(opt.sGmax)
        opt.sGmax = sg;
    end
    
    drho = fluid.rhoW(p) - fluid.rhoG(p);
    sg = free_sg(sg, opt.sGmax, fluid.res_water, fluid.res_gas);

    SP = compute_untruncated_SP(H, drho, sg, table) ./ opt.kscale;
    
    pc = interpTable(table.SP, table.p, SP) .* opt.kscale;
    
end

% ----------------------------------------------------------------------------
function [SP, H_trunc] = compute_untruncated_SP(H, drho, sg, table)


    drho_g = drho * norm(gravity);
    dP_aquifer = H .* drho_g;
    
    SP = sg .* dP_aquifer;
    
    % compute the full height of the plume represented by SP
    p = interpTable(table.SP, table.p, SP);
    H_fullplume = p ./ drho_g;
    
    x = max(H_fullplume - H, 0); % how much the full plume extends below aquifer bottom
    trunc_ind = x > 0;
    H_trunc = 0 * H;
    if isa(sg, 'ADI')
        H_trunc = double2ADI(H_trunc, sg);
    end
    
    if any(trunc_ind)
        %disp('Truncated columns present'); % @@ This line can be commented out.
        % there were columns with truncated plumes.  For these columns, compute the
        % value of SP that corresponds to a full, untruncated column.
        
        % determine local variables that only relate to truncated columns
        x_loc = x(trunc_ind);
        drho_g_loc = drho_g(trunc_ind);
        SP_loc = SP(trunc_ind);
        H_loc = H(trunc_ind);

        % determine capillary pressure, 'pb', at aquifer bottom
        f = @(pb) interpTable(table.p, table.SP,  pb) + SP_loc - ...
                  interpTable(table.p, table.SP, pb + drho_g_loc .* H_loc);
        df = @(pb, ixs) dinterpTable(table.p, table.SP,  value(pb(ixs))) - ...
                   dinterpTable(table.p, table.SP, value(pb(ixs) + drho_g_loc(ixs) .* H_loc(ixs))); 
        
        MAX_ITER = 2000; % should be largely enough 
        TOL = 1e-3 * max(dP_aquifer);

        pb = x_loc .* drho_g_loc; % initial guess
        
        for iter = 1:MAX_ITER
            res = f(pb);
            if any(res > TOL)
                ixs = res > TOL;
                pb(ixs) = pb(ixs) - res(ixs) ./ df(pb, ixs); % apply Newton
            else
                % all bottom capillary pressures are within tolerance
                break;
            end
        end
        if any(res > TOL)
            warning('Did not converge when computing untruncated SP');
        end

        SP(trunc_ind) = SP_loc + interpTable(table.p, table.SP, pb);
        H_trunc(trunc_ind) = pb ./ drho_g_loc;
        
    end
end

% ----------------------------------------------------------------------------
function pc = pcWG_htable(sg, table, fluid, H, varargin)
% note that this function is strictly valid only for constant H and incompressible fluid

    opt = merge_options(struct('sGmax', []), varargin{:});
    if isempty(opt.sGmax)
        opt.sGmax = sg;
    end

    SH = free_sg(sg, opt.sGmax, fluid.res_water, fluid.res_gas) .* H;
    h = interpTable(table.SH, table.h, SH);

    drho = fluid.rhoWS - fluid.rhoGS;
    pc = h .* drho * norm(gravity);
end

% ----------------------------------------------------------------------------
function kr = krG_htable(sg, table, fluid, H, varargin)
    opt = merge_options(struct('sGmax', []), varargin{:});
    if isempty(opt.sGmax)
        opt.sGmax = sg;
    end
    
    SH = free_sg(sg, opt.sGmax, fluid.res_water, fluid.res_gas) .* H;
    assert(all(SH ./ H <= 1));
    
    gh = interpTable(table.SH, table.h, SH); % h corresponding to our SH
    kr = interpTable(table.h, table.krH, gh) / H; % averaged relperm corresponding to our SH
end


% ----------------------------------------------------------------------------
function issue_warnings(type, fluid, Gt, rock)
    
    switch type
      case 'S table'
        if max(Gt.cells.H) > min(Gt.cells.H)
            warning(['Provided grid breaks chosen relperm model (S table) ' ...
                     'assumption of constant aquifer thickness.']);
        end
        if fluid.rhoG(1e5) ~= fluid.rhoG(1e6) % test if fluid compresses
            warning(['Provided fluid breaks chosen relperm model (S table) ' ...
                     'assumption of incompressible fluid.']);
        end
      case 'P-scaled table'
        if max(rock.poro) ~= min(rock.poro)
            warning(['The ''P-scaled table'' relperm model is used, which ' ...
                     'assumes the same capillary pressure function ' ...
                     'everywhere.  However, provided rock properties are ' ...
                     'heterogeneous.  Consider using the k-scaled model.']);
        end
      case {'P-K-scaled table', 'linear cap.'}
        % nothing 
      otherwise
        error('Unknown VE relperm model.');
    end
end

% ----------------------------------------------------------------------------
function Pmax = determine_upper_cap_press_limit(fluid, Gt, invPc3D)
    dp = 1 * Pascal;
    max_dp = 2 * mega * Pascal; % we do not foresee upscaled capillary pressures
                                % reaching higher than this in practice.
    s = invPc3D(dp);
    tol = 1e-3;
    while (s < max_dp) && (invPc3D(dp) > fluid.res_water + tol)
        dp = dp * 2;
    end
    
    % Density difference at surface (or reference) level. Expected to generally
    % decrease in depth, so this serves as an upper estimate of density
    % difference.
    drho = fluid.rhoWS - fluid.rhoGS; 
                                      
    % maximum drop in capillary pressure along vertical direction in aquifer
    max_cap_drop = max(Gt.cells.H) * drho * norm(gravity);
    
    % ensures that a cap. pressure of 'Pmax' at caprock level is enough to yield a
    % fully saturated vertical column of CO2 even at the thickest point of the
    % aquifer
    Pmax = dp + max_cap_drop; 
end

% ----------------------------------------------------------------------------
function table = make_CO2_table_h_based(invPc3D, kr3D, Gt, samples, Pmax, drho_surf)
    
    Hmax = Pmax / drho_surf / norm(gravity); % NB: using h directly is only strictly
                                             % valid in the incompresssible case
    h = linspace(0, Hmax, samples);
    dh = h(2) - h(1);
    
    % fine scale saturation and relperm at vertical % distance 'h' from the plume tip
    sw_h = invPc3D(h * drho_surf * norm(gravity));  % water saturation (fine-scale)
    sg_h = 1 - sw_h;                                 % CO2 saturation (fine-scale)
    kr_h = kr3D(sg_h);
    
    % vertical integrated values
    SH  = (cumsum(sg_h) - sg_h / 2) * dh;
    krH = (cumsum(kr_h) - kr_h / 2) * dh;
    
    table.SH  = truncate_at_aquifer_bottom(SH, h, max(Gt.cells.H));
    table.krH = truncate_at_aquifer_bottom(krH, h, max(Gt.cells.H));
    table.h = h;
end

% ----------------------------------------------------------------------------
function table = make_CO2_table_p_based(invPc3D, kr3D, samples, Pmax)
    
    p = linspace(0, Pmax, samples);
    dp = p(2) - p(1);
    
    % fine scale saturation and relperm at vertical position with capillary pressure 'p'
    sw_p = invPc3D(p);
    sg_p = 1 - sw_p;
    kr_p = kr3D(sg_p);
    
    % With nonzero entry pressure, more than one of the first saturation
    % values will likely be zero.  To ensure an invertible function, we
    % smooth this part.
    last_zero_s = find(sg_p == 0, 1, 'last');
    SMOOTH_FAC = 0.1;
    if last_zero_s > 1
        sg_p(1:last_zero_s) = linspace(0, SMOOTH_FAC * sg_p(last_zero_s + 1), last_zero_s);
    end
    
    % vertical integrated values
    table.SP  = (cumsum(sg_p) - sg_p / 2) * dp;
    table.krP = (cumsum(kr_p) - kr_p / 2) * dp;
    table.p   = p;
    
end

% ----------------------------------------------------------------------------
function vec_truncated = truncate_at_aquifer_bottom(vec, h, H)
    
    truncate = h > H;
    
    % 'chop off' the (theoretical) part of the plume that extends below aquifer
    % bottom
    vec_truncated = vec;
    vec_truncated(trunctate) = vec(truncate) - interpTable(h, vec, h(truncate) - H);
end
