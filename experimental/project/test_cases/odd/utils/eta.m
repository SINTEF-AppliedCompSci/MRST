function [val, int_val] = eta(Gcomp, EOS, P, T, Tgrad, dz, order)
%
% SYNOPSIS:
%   function [val, int_val] = eta(P, T, Tgrad, dz, order)
%n
% PARAMETERS:
%   Gcomp - gravity component along z              (typically 'gravity * cos (theta))
%   EOS   - Equation of state                      (typically from CO2props)
%   P     - Pressure at reference level            (vector of values; can be ADI)
%   T     - Temperature at reference level         (vector of values; can be ADI)
%   Tgrad - Temperature gradient                   (scalar)
%   dz    - Extrapolation distance in z -direction (scalar)
%   order - order of extrapolation                 (1: linear, 2: quadratic)          
%
% RETURNS:
%   val     - estimated value of 'eta' (estimated ratio of density value at
%             vertical distance dz from reference point, over reference point
%   int_val - dz-averaged value of 'eta, i.e. the integral of 'eta from 0 to
%             dz, divided by dz
%
% @@TODO: 
%   Add support for jump across liquid-vapor boundary!
%

if order == 0
    val     = 0*P;     % This roundabout way of setting 'val' to a vector of
    val     = val + 1; % ones is to ensure compatibility with ADI framework. 
    int_val = val; 
    return
end

% -- If we got here, approximation should at least be linear. -- 

% computing reference density
rhoRef  = EOS.rho(P, T);
val     = rhoRef;
int_val = val;

% adding linear term
dP       = Gcomp * rhoRef;
dT       = Tgrad;
rho_dP   = EOS.rhoDP(P, T);
rho_dT   = EOS.rhoDT(P, T);
d_dz_rho = (dP .* rho_dP + dT .* rho_dT);

val      = val + (d_dz_rho .* dz);
int_val  = int_val + (1/2) * (d_dz_rho .* dz);

if (order == 1)
    val     = val ./rhoRef;
    int_val = int_val ./rhoRef;
    return
end

% -- If we got here, approximation should be quadratic --

assert(order == 2); % other orders not supported
    
% Adding quadratic term
dPP = Gcomp .* d_dz_rho;

d2_dz2_rho = (dP.^2 .* EOS.rhoDPP(P, T)            + ...
                 2 .* dP .* dT .* EOS.rhoDPT(P, T) + ...
                 dT.^2 .*EOS.rhoDTT(P, T)          + ...
                 dPP .* rho_dP);

val     = (val + (0.5 * d2_dz2_rho .* dz .* dz)) ./ rhoRef;
int_val = (int_val + (1/6) * d2_dz2_rho .* dz .* dz) ./ rhoRef;