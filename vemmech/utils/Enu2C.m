function C = Enu2C(E, nus, G)
%
%
% SYNOPSIS:
%   function C = Enu2C(E, nus, G)
%
% DESCRIPTION: For each cell, construct the 3x3 (in 2D) or 6x6 (in 3D) matrix
% representing the elasticity tensor.
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
% PARAMETERS:
%   E   - Young's modulus (one entry per cell)
%   nus - Poisson ratio (one entry per cell)
%   G   - Grid
%
% RETURNS:
%   C - (k,n) matrix, where k=3^2 (2D) or k=6^2 (3D), and n is the number of
%       cells.  Each row thus represents the entries of the elasticity tensor
%       for a specific cell.

% convert from Voigt notation to matrix in VEM format of vem formulation

   if(G.griddim == 2)
      z = zeros(numel(nus), 1); 
      C = [reshape([1 - nus, nus    , z          ]', [], 1), ...
           reshape([nus    , 1 - nus, z          ]', [], 1), ...
           reshape([z      , z      , 1 - 2 * nus]', [], 1) / 2]; 
      nlin = 3; 
   else
      assert(G.griddim == 3); 
      nlin = 6; 
      zzz = zeros(numel(nus), 3); 
      zz = zeros(numel(nus), 2); 
      z = zeros(numel(nus), 1); 
      C = [reshape([1 - nus, nus        , nus    , zzz  ]', [], 1)    , ...
           reshape([nus    , 1 - nus    , nus    , zzz  ]', [], 1)    , ...
           reshape([nus    , nus        , 1 - nus, zzz  ]', [], 1)    , ...
           reshape([zzz    , 1 - 2 * nus, zz            ]', [], 1) / 2, ...
           reshape([zzz    , z          , 1 - 2 * nus, z]', [], 1) / 2, ...
           reshape([zzz    , zz         , 1 - 2 * nus   ]', [], 1) / 2]; 
      
   end
   fac = (E ./ ((1 + nus) .* (1 - 2 * nus))); 
   C   = reshape(C', nlin * nlin, [])';
   C   = bsxfun(@times, C, fac);
end

