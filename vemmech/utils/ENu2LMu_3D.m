function [lambda, mu] = ENu2LMu_3D(E, nu)
%
%
% SYNOPSIS:
%   function [lambda, mu] = ENu2LMu_3D(E, nu)
%
% DESCRIPTION: Compute the Lamé parameters lambda and mu in the 2D case from the Young's
% modulus E and Poisson's ratio Nu. In the orthogonal direction (with respect
% to the 2D plane), it is assumed zero displacement. In particular, it
% implies that Nu = 0.5 corresponds as in the 3D case to an incompressible
% material. An other assumption, which is not used here, could be that there
% is no stress in the orthogonal direction.
%
% PARAMETERS:
%   E  - Young's modulus
%   nu - Poisson's ratio
%
% RETURNS:
%   lambda - First Lamé coefficient
%   mu     - Second Lamé coefficient
%
% EXAMPLE:
%
% SEE ALSO:
%

   lambda = (E .* nu) ./ ( (1 + nu) .* (1 - 2 * nu));
   mu = E ./ (2 * (1 + nu));

end
