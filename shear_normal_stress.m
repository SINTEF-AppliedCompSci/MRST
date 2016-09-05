function constit = shear_normal_stress(Nc, Nd, mu, lambda, phi)
% Define stiffness matrix for linear elastic medium.
%
% Parameters
% Nc: Number of cells in the grid
% Nd: Number of dimensions
% mu: First lame parameter, cellwise (Nc x 1 vector)
% lambda: Second Lame parameter, cellwise (Nc x 1 vector)
% phi: Extra parameter, never been used, stay away unless you know what you
%   are doing 
% 
% Returns:
%   Constitutive relation, one matlab cell per computational cell.
%
%{
Copyright 2015-2016, University of Bergen.

This file is part of FVBiot.

FVBiot is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

FVBiot is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this file.  If not, see <http://www.gnu.org/licenses/>.
%} 

alpha = 0; %1e-10; 
constit = cell(Nc,1);

tmpInd = zeros(Nd,Nd);
for iter1 = 1 : Nd
    for iter2 = 1: Nd
        tmpInd(iter1,iter2) = sub2ind_loc([Nd,Nd],iter1, iter2);
    end
end

for iter1 = 1:Nc;
    constit{iter1} = zeros(Nd^2);
    for iter3 = 1:Nd %iter3 indicates the normal vector of the stress surface
        for iter4 = 1:Nd %iter4 indicates the unit directions of stress
            %The constitutive matrix gives the directional stress due to a given
            %displacement gradient
            if (iter3==iter4)
                temp = lambda(iter1)*eye(Nd,Nd) + phi(iter1)*(ones(Nd)-eye(Nd));
            else
                temp = phi(iter1)*eye(Nd);
            end
            if 1 % Nonsymmetric formulation
                temp(iter3,iter4)=temp(iter3,iter4)+(1+alpha)*mu(iter1);
                temp(iter4,iter3)=temp(iter4,iter3)+(1-alpha)*mu(iter1);
            else % Imposes pointwise symmetry - makes MPFA local system underdetermined:
                temp(iter3,iter4)=temp(iter3,iter4)+mu(iter1);
                temp(iter4,iter3)=temp(iter4,iter3)+mu(iter1);
            end
            constit{iter1}(tmpInd(iter3, iter4),:) = reshape(temp,1,Nd^2);
        end
    end
end
