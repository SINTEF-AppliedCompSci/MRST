function [val avg_val] = fp_eta(Gcomp, EOS, P, T, Tgrad, dz, order)%
%
% Compute the product of the functions 'eta' and 'fp' for a given depth 'dz'. 
% 'eta' is defined as the ratio of density at 'dz' over top density, while
% 'fp' is defined as the derivative of pressure at 'dz' as a function of top
% pressure.  This product is useful in the computation of mass flux as depth
% 'dz'. 
%
% SYNOPSIS:
%   function [val avg_val] = fp_eta(Gcomp, EOS, P, T, Tgrad, dz, order)
%
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
%   val     - estimation of product of 'eta' and 'fp' at depth 'dz'
%   avg_val - estimation of the dz-averaged values of the product of 'eta'
%             and 'fp', from the top down to 'dz'. 
%
% EXAMPLE:
%
% SEE ALSO:
%

%% Computing zeroth term

eta     = 0*P;     % this roundabout way of setting 'val' to a vector of ones
eta     = eta+1;   % is to ensure compatibility with the ADI framework       
fp      = eta;
val     = eta.*fp; % ones x ones = ones (but written this way for clarity)
avg_val = val;     % val is constant, so it is equal to its average

if order == 0
    return;
end

%% -- If we got here, approximation should at least be linear. --

%% computing reference density
rhoRef  = EOS.rho(P, T);

%% adding linear term

% computing d_dz_eta
dP       = Gcomp * rhoRef;
dT       = Tgrad;
rho_dP   = EOS.rhoDP(P, T);
rho_dT   = EOS.rhoDT(P, T);
d_dz_eta = (dP .* rho_dP + dT .* rho_dT)./rhoRef;

% computing d_dz_fp
beta = EOS.beta(P, T);
d_dz_fp = dP .* beta;

% domputing d_dz_eta_fp
d_dz_eta_fp = d_dz_fp + d_dz_eta; % remember that reference fp = eta = 1

% updating fp_eta and its average with linear term
val = val + d_dz_eta_fp .* dz; 
avg_val = avg_val + (1/2) * d_dz_eta_fp .* dz;

if order == 1
    return;
end

%% -- If we got here, approxmation should be quadratic. --
assert(order == 2); % higher orders not supported

%% adding quadratic term

% computing d2_dz2_eta
d_dz_rho = (dP .* rho_dP + dT .* rho_dT);
dPP = Gcomp * d_dz_rho;
d2_dz2_eta = (dP.^2 .* EOS.rhoDPP(P, T)            + ...
                 2 .* dP .* dT .* EOS.rhoDPT(P, T) + ...
                 dT.^2 .*EOS.rhoDTT(P, T)          + ...
                 dPP .* rho_dP) ./ rhoRef;

% computing d2_dz2_fp
beta_Two = EOS.rhoDPP(P, T) ./ rhoRef;
chi      = EOS.rhoDPT(P, T) ./ rhoRef;
d2_dz2_fp = dP .* (dP .* (beta.^2 + beta_Two) + dT .* chi);

% computing d2_dz2_eta_fp
d2_dz2_eta_fp = d2_dz2_eta + (2 * d_dz_eta .* d_dz_fp) + d2_dz2_fp;

% updating fp_eta and its average with quadratic term
val = val + (1/2) * d2_dz2_eta_fp .* (dz.^2);
avg_val = avg_val + (1/6) * d2_dz2_eta_fp .* (dz.^2);
