function p = hydrostatic(Gt, cell_ref, at_top, rho_ref, beta_ref, p_ref, theta, slopedir)
% +++++++DEPRECATED+++++++ (WE USE CONSTANT DENSITY FOR WATER INSTEAD, ATM)
% Compute a hydrostatic pressure field over the vertically-averaged grid.
%
% SYNOPSIS:
%   function p = hydrostatic(Gt, cell_ref, rho_ref, beta_ref, theta, slopedir, at_top)
%
% PARAMETERS:
%   Gt       - Vertically averaged grid
%   cell_ref - index to cell for which reference values are provided
%   rho_ref  - density (of water) at reference cell
%   beta_ref - compressibility (of water) at reference cell
%   p_ref    - pressure (of water) at reference cell
%   theta    - dip angle of aquifer
%   slopedir - dip direction
%   at_top   - 'true' if hydrostatic pressure is to be computed for the top
%              (otherwise, it will be for the bottom).
%
% RETURNS:
%   p - hydrostatic pressure field over the cells of Gt
%

% computing height differential between each cell in Gt and reference cell,
% disregarding dip angle (NB: remember inverted z-axis)

if (at_top)
    ref_z = Gt.cells.z(cell_ref);
    dz = Gt.cells.z - ref_z;
else
    ref_z = Gt.cells.z(cell_ref) + Gt.cells.H(cell_ref);
    dz = (Gt.cells.z + Gt.cells.H) - ref_z;
end

% including dip angle in height differential
slopedir = slopedir ./ norm(slopedir);   % ensuring normalized vector
ref_dist = Gt.cells.centroids - Gt.cells.centroids(cell_ref); % list of 2D vectors
dz = dz - sin(theta) * ref_dist * slopedir'; % moving in the direction of
                                             % 'slopedir' takes us to smaller
                                             % depths (i.e. smaller values of 'z')

% computing 'linear pressure differential'
rgh = cos(theta) * norm(gravity) * rho_ref * dz;

% computing hydrostatic pressure for each cell analytically (valid under the
% hypothesis of constant compressibility)
if beta_ref > 0
    % compressible water
    p = p_ref - (1/beta_ref) * log(1 - beta_ref * rgh);
else
    % incompressible water
    p = p_ref + rgh;
end
