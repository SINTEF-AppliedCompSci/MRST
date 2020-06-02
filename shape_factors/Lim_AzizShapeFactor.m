classdef Lim_AzizShapeFactor < ShapeFactor
%
% SYNOPSIS:
%   model = SimpleTransferFunction(shape_factor_name, block_dimension)
%
% DESCRIPTION: 
%   Single phase transfer function as originally described in Warren and 
%   Root (1963)
%
% PARAMETERS:
%   shape_factor_name  - Name of shape factore (see dual-porosity module)
%   block_dimension    - Dimensions of dual-continuum block (defined by 
%                        fracture spacing)
%
% RETURNS:
%   class instance
%
% EXAMPLE:
%
% SEE ALSO: TransferFunction in dual-porosity
%
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

    properties

    end

    methods
        function shapefactor = Lim_AzizShapeFactor(block_dimension)
          shapefactor = shapefactor@ShapeFactor(block_dimension);
        end

        function [sigma] = calculate_shape_factor(sf,model)

            kx = model.rock_matrix.perm(:,1);
            try
                ky = model.rock_matrix.perm(:,2);
                kz = model.rock_matrix.perm(:,3);
            catch
                ky = model.rock_matrix.perm(:,1);
                kz = model.rock_matrix.perm(:,1);
            end
            
            if size(sf.block_dimension, 2) < 3
                lx = sf.block_dimension(:,1);
                ly = sf.block_dimension(:,2);
                sigma=((pi^2)./(kx.*ky).^(1/2))...
                      .*(kx./(lx.^2)+ky./(ly.^2));
            else
                lx = sf.block_dimension(:,1);
                ly = sf.block_dimension(:,2);
                lz = sf.block_dimension(:,3);
                sigma=((pi^2)./(kx.*ky.*kz).^(1/3))...
                      .*(kx./(lx.^2)+ky./(ly.^2)+kz./(lz.^2));
            end
        end
    end
end
