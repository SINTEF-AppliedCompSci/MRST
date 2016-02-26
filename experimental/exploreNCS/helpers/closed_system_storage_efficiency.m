function E = closed_system_storage_efficiency( cw, cr, P_over, P_init )
% Compute storage efficiency factor for a closed-system
% (compressibility-limited)

% This equation comes from Zhou et al. 2008: "A method for quick
% assessment of CO2 storage capacity in closed and semi-closed saline
% formations"

% cw is water compressibility, Pa^-1

% cr is rock (pore volume) compressibility
%   - default is 1e-5 / barsa (gets converted to Pa^-1), which is the
%   pvMult_fac used in makeVEFluid.

assert( all(P_over > P_init) );
E = (cw + cr) * min(P_over - P_init); % unitless
%E = (cw + cr) * mean(repmat(min(P_over), numel(P_init), 1) - P_init); % unitless
%E = (cw + cr) * mean(P_over - P_init); % unitless


% NB: the minimum P_over is taken in order to get a conservative
% estimate for E


end

