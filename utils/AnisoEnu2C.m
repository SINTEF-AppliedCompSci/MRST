function [C, invC, invCi] = AnisoEnu2C(E, nu, mu, G, varargin)
%
% SYNOPSIS:
%   function [C, invC, invCi] = AnisoEnu2C(E, nu, mu, G)
%
% DESCRIPTION: 
%   Function to compute cellwise, 3x3 (in 2D) or 6x6 (in 3D), anisotropic 
%   stiffness tensors (limited to ORTHOTROPIC materials and subsets
%   therein). Compliance tensors are also computed. 
%
%   Note: Notation used herein follows notation on orthotropic tensors 
%   according to Vannucci (2018). Anisotropic elasticity. Singapore: Springer.
%
% PARAMETERS:
%   G   - Grid struc
%   C_m - matrix stiffness tensor
%   C_f - fracture stiffness tensor
%   varargin - stiffness tensor calculation options
%
% RETURNS:
%   C (stiffness), invC (compliance), invCi (trace of compliance)
%
% EXAMPLE:
%
% SEE ALSO:
%   Enu2C
%{
Copyright 2009-2020 SINTEF ICT, Applied Mathematics.

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

% For each cell, construct the 3x3 (in 2D) or 6x6 (in 3D) matrix representing
% the anisotropic elasticity tensor (limited to ORTHOTROPIC materials and subsets
% therein)
%
% SYNOPSIS:
%   function C = AnisoEnu2C(E, nu, mu, G)
%
% Notation used herein follows notation on orthotropic tensors according to
% Vannucci, P., 2018. Anisotropic elasticity. Singapore: Springer.

% 2D first
   if(G.griddim == 2)      
       [mu12] = deal(mu(:,1));
       z = zeros(numel(mu12), 1); 
       nlin = 3;
       id = [1; 1; 0];
      
       if nargin > 4          
           if strcmp(varargin{1}, 'plane_stress')
               [E1,E2] = deal(E(:,1),E(:,2));
               [nu21] = deal(nu(:,1));
               nu12 = (nu21.*E1)./E2;
               C = [reshape([        E1./(1-nu12.*nu21), (nu21.*E1)./(1-nu12.*nu21),        z]', [], 1), ...
                    reshape([(nu12.*E2)./(1-nu12.*nu21),         E2./(1-nu12.*nu21),        z]', [], 1), ...
                    reshape([                         z,                          z,     mu12]', [], 1)];

               invC = [reshape([    1./E1,  -nu21./E2,            z]', [], 1), ...
                       reshape([-nu12./E1,      1./E2,            z]', [], 1), ...
                       reshape([        z,          z,      1./mu12]', [], 1) ];

               invCi = invC*id;
               invC  = reshape(invC', nlin*nlin, [])';
               invCi = reshape(invCi, nlin , [])';

           elseif strcmp(varargin{1}, 'plane_strain')
               [E1,E2,E3] = deal(E(:,1), E(:,2), E(:,3));
               [nu21,nu31,nu32] = deal(nu(:,1),nu(:,2),nu(:,3));
               nu12 = (nu21.*E1)./E2;
               nu13 = (nu31.*E1)./E3;
               nu23 = (nu32.*E2)./E3;
               del = (1 - nu12.*nu21 - nu23.*nu32 - nu13.*nu31 - 2.*nu21.*nu32.*nu13);
 
               C = [reshape([   ((1 - nu23.*nu32)./del).*E1,  ((nu21 + nu23.*nu31)./del).*E1,        z]', [], 1), ...
                    reshape([((nu12 + nu32.*nu13)./del).*E2,       ((1-nu13.*nu31)./del).*E2,        z]', [], 1), ...
                    reshape([                             z,                               z,     mu12]', [], 1)];

               invC = [reshape([                  (1-nu13.*nu31)./E1,   -nu21.*(nu12+nu32.*nu13)./(nu12.*E2),            z]', [], 1), ...
                       reshape([-nu12.*(nu21+nu23.*nu31)./(nu21.*E1),                      (1-nu23*nu32)./E2,            z]', [], 1), ...
                       reshape([                                   z,                                      z,      1./mu12]', [], 1) ];

               invCi = invC*id;
               invC  = reshape(invC', nlin*nlin, [])';
               invCi = reshape(invCi, nlin , [])';     
           end
       else
       % assume plane stress    
           [E1,E2] = deal(E(:,1),E(:,2));
           [nu21] = deal(nu(:,1));
           nu12 = (nu21.*E1)./E2;
           C = [reshape([        E1./(1-nu12.*nu21), (nu21.*E1)./(1-nu12.*nu21),        z]', [], 1), ...
                reshape([(nu12.*E2)./(1-nu12.*nu21),         E2./(1-nu12.*nu21),        z]', [], 1), ...
                reshape([                         z,                          z,     mu12]', [], 1)];
   
           invC = [reshape([    1./E1,  -nu21./E2,            z]', [], 1), ...
                   reshape([-nu12./E1,      1./E2,            z]', [], 1), ...
                   reshape([        z,          z,      1./mu12]', [], 1) ];

           invCi = invC*id;
           invC  = reshape(invC', nlin*nlin, [])';
           invCi = reshape(invCi, nlin , [])';
       end
   % 3D 
   else
       assert(G.griddim == 3);
       [E1, E2, E3] = deal(E(:,1), E(:,2), E(:,3));
       [nu21, nu31, nu32] = deal(nu(:,1) ,nu(:,2), nu(:,3));
       [mu12, mu13, mu23] = deal(mu(:,1), mu(:,2), mu(:,3));
      
       % Calculate missing parameters by taking advantage of the symmetry
       nu12 = (nu21.*E1)./E2;
       nu13 = (nu31.*E1)./E3;
       nu23 = (nu32.*E2)./E3;
       zzz = zeros(numel(nu21), 3);
       zz = zeros(numel(nu21), 2);
       z = zeros(numel(nu21), 1);
       nlin = 6;
       id = [1; 1; 1; 0; 0; 0];
      
       % Caluclate elasticity and compliance tensors
       del = (1 - nu12.*nu21 - nu23.*nu32 - nu13.*nu31 - 2.*nu21.*nu32.*nu13);
      
       C = [reshape([   ((1 - nu23.*nu32)./del).*E1,  ((nu21 + nu23.*nu31)./del).*E2,    ((nu31 + nu21.*nu32)./del).*E3,                    zzz]', [], 1) , ...
            reshape([((nu12 + nu32.*nu13)./del).*E2,       ((1-nu13.*nu31)./del).*E2,    ((nu32 + nu12.*nu31)./del).*E3,                    zzz]', [], 1) , ...
            reshape([((nu13 + nu12.*nu23)./del).*E3,  ((nu23 + nu21.*nu13)./del).*E3,         ((1-nu12.*nu21)./del).*E3,                    zzz]', [], 1) , ...
            reshape([                                                                                                           zzz,   mu23,             zz]', [], 1) , ...
            reshape([                                                                                                           zzz,      z,   mu13,      z]', [], 1) , ...
            reshape([                                                                                                           zzz,             zz,   mu12]', [], 1) ];

       invC = [reshape([1./E1    ,   -nu21./E2,  -nu31./E3,                           zzz]', [], 1), ...
               reshape([-nu12./E1,       1./E2,  -nu32./E3,                           zzz]', [], 1), ...
               reshape([-nu13./E1,   -nu23./E2,      1./E3,                           zzz]', [], 1), ...
               reshape([                               zzz,  1./mu23,                  zz]', [], 1), ...
               reshape([                               zzz,        z,  1./mu13,         z]', [], 1), ...
               reshape([                               zzz,                 zz,   1./mu12]', [], 1)];

       invCi = invC*id;
       invCi = reshape(invCi, nlin , [])';
       invC  = reshape(invC', nlin*nlin, [])';
   end

   C   = reshape(C', nlin * nlin, [])';
end
