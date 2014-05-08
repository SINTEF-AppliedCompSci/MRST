function pdiff = pressureDifference(fluid, dh, P, T, theta, P_is_ref_pressure)
    
% Helper function used by both Pt2PbPi and Pb2PtPi
    
% NB: In the case of fluids with 'HORIZONTAL' compressibility, the density does
% not vary in the vertical direction.  However, in order to choose the
% density to use for a vertical column, a reference pressure must be used.
% This can either be 'P', or 'P' + 'pdiff'.  If 'P' is the reference
% pressure, the last argument ('P_is_ref_pressure') should be set to 'true',
% otherwise, it should be set to false.  If the fluid has 'FULL'
% compressibility, or is 'INCOMPRESSIBLE', this last argument does not matter.
    
    comp = fluid.compressible;
    gval = norm(gravity) * cos(theta);    
    
    % Check if we should do a correction in pressure to account for vertical
    % density changes
    if strcmp(comp,'FULL')
        poly = VEpolys;

        beta0 = fluid.beta(P, T);
        bder0 = fluid.bder(P, T);
        rho0  = fluid.rho(P, T);
        
        pdiff = rho0 .* gval .* poly.intAlpha(dh, rho0, beta0, bder0, theta);
        
    elseif strcmp(comp, 'INCOMPRESSIBLE')
        
        % no vertical correction in pressure, and reference density is known (constant)
        rho0 = fluid.rho(P, T); % just dummy arguments here!
        pdiff = rho0 .* gval .* dh;
        
    else
        assert (strcmp(comp, 'HORIZONTAL'));
        
        % no vertical correction in pressure (constant rho).  Do we know the
        % reference pressure/rho?
        if P_is_ref_pressure
            % We know the reference pressure.  The rest is easy.
            rho0 = fluid.rho(P, T);
            pdiff = rho0 .* gval .* dh;
        else
            % We do not know the reference pressure (P_ref).  We have to use
            % a root-finding approach to reconstruct it, based on the
            % relation: P_ref = P + rho(P_ref) * gval * dh
            pref = arrayfun(@(p, h) compute_ref_pressure(p, h, fluid, T, gval), P, dh);
            pdiff = pref-P;  % we need to return the pressure _difference_
        end
    end
end
%%                                                                              

function pref = compute_ref_pressure(p, h, fluid, T, gval)
%  Function to compute 'pref' such that:
%  pref = P + fluid.rho(pref, T) * gval * h
    pref = fzero(@(x) residual(x, p, h, fluid, T, gval), p); 
end
%%                                                                              

function res = residual(x, p, h, fluid, t, gval)
    res = (x - p) - fluid.rho(x, t) * gval * h;
end
%%                                                                              

