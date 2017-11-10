function [C, invC, invCi] = Enu2C(E, nu, G)
%
% For each cell, construct the 3x3 (in 2D) or 6x6 (in 3D) matrix representing
% the elasticity tensor.
%
% SYNOPSIS:
%   function C = Enu2C(E, nu, G)
%
% DESCRIPTION:
%
% In 3D, the matrix of the elasticity tensor for a given cell is written:
%
% |  1-nu  nu    nu    0          0          0         |
% |  nu    1-nu  nu    0          0          0         |
% |  nu    nu    1-nu  0          0          0         |            E
% |  0     0     0     (1-2nu)/2  0          0         |  x  --------------
% |  0     0     0     0          (1-2nu)/2  0         |     (1+nu) (1-2nu)
% |  0     0     0     0          0          (1-2nu)/2 |
%
%
% In 3D, the inverse of the elasticity tensor is given by
%  [ [   1 / E,   -(nu / E), -(nu / E),      0,            0,            0       ]
%    [ -(nu / E),   1 / E,   -(nu / E),      0,            0,            0       ]
%    [ -(nu / E), -(nu / E),   1 / E,        0,            0,            0       ]
%    [     0,         0,         0,     2*(nu + 1) / E,      0,            0       ]
%    [     0,         0,         0,          0,       2*(nu + 1) / E,      0       ]
%    [     0,         0,         0,          0,            0,       2*(nu + 1) / E ] ]
%
% In 2D, the inverse is given
%    [ [ 1 - nu,  -nu,   0 ]            1 + nu
%      [  -nu,   1 - nu, 0 ]     x  --------------
%      [   0,      0,    2 ] ]            E
%
% The factors 2 in the expressions of the inverse above are due to Voigts
% notations, see https://en.wikipedia.org/wiki/Voigt_notation
%
%
% invCi = C^{-1}*(Identity tensor)
%
% PARAMETERS:
%   E   - Young's modulus (one entry per cell)
%   nu - Poisson ratio (one entry per cell)
%   G   - Grid
%
% RETURNS:
%   C - (k,n) matrix, where k=3^2 (2D) or k=6^2 (3D), and n is the number of
%       cells.  Each row thus represents the entries of the elasticity tensor
%       for a specific cell.


   if(G.griddim == 2)
      z = zeros(numel(nu), 1);
      o = ones(numel(nu), 1);
      C = [reshape([1 - nu, nu    ,        z]', [], 1), ...
           reshape([nu    , 1 - nu,        z]', [], 1), ...
           reshape([      z,     z, 1 - 2*nu]', [], 1) / 2];

      invC = [reshape([1 - nu,    -nu,   z]', [], 1), ...
              reshape([   -nu, 1 - nu,   z]', [], 1), ...
              reshape([     z,      z, 2*o]', [], 1) ];

      id = [1; 1; 0];

      invCi = invC*id;

      nlin = 3;

      invfac = (1 + nu)./E;
      invC  = reshape(invC', nlin*nlin, [])';
      invC = bsxfun(@times, invC, invfac);

      invCi = reshape(invCi, nlin , [])';
      invCi = bsxfun(@times, invCi, invfac);

   else
      assert(G.griddim == 3);
      nlin = 6;
      zzz = zeros(numel(nu), 3);
      zz = zeros(numel(nu), 2);
      z = zeros(numel(nu), 1);
      C = [reshape([1 - nu, nu    , nu    ,                                zzz]', [], 1)    , ...
           reshape([nu    , 1 - nu, nu    ,                                zzz]', [], 1)    , ...
           reshape([nu    , nu    , 1 - nu,                                zzz]', [], 1)    , ...
           reshape([                   zzz, 1 - 2 * nu,                     zz]', [], 1) / 2, ...
           reshape([                   zzz,          z, 1 - 2 * nu,          z]', [], 1) / 2, ...
           reshape([                   zzz,                     zz, 1 - 2 * nu]', [], 1) / 2];

      invC = [reshape([1./E  , -nu./E, -nu./E,                                         zzz]', [], 1), ...
              reshape([-nu./E, 1./E  , -nu./E,                                         zzz]', [], 1), ...
              reshape([-nu./E, -nu./E, 1./E  ,                                         zzz]', [], 1), ...
              reshape([                   zzz, 2*(1 + nu)./E,                           zz]', [], 1), ...
              reshape([                   zzz,             z, 2*(1 + nu)./E,             z]', [], 1), ...
              reshape([                   zzz,                           zz, 2*(1 + nu)./E]', [], 1)];

      id = [1; 1; 1; 0; 0; 0];

      invCi = invC*id;
      invCi = reshape(invCi, nlin , [])';
      invC  = reshape(invC', nlin*nlin, [])';

   end
   fac = (E ./ ((1 + nu) .* (1 - 2 * nu)));
   C   = reshape(C', nlin * nlin, [])';
   C   = bsxfun(@times, C, fac);
end
