classdef SimpleTransferFunction < TransferFunction
%
% SYNOPSIS:
%   transferfunction = SimpleTransferFunction(shape_factor_name, block_dimension)
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
        shape_factor_object
        shape_factor_name
    end

    methods

        function transferfunction = SimpleTransferFunction(shape_factor_name, block_dimension)

            transferfunction = transferfunction@TransferFunction();
            transferfunction.nphases = 1; % Although this is a single phase transfer function,  
                                          % parent class requires us to set nphases = 2 to 
                                          % check for correct fields
                                          
            if (nargin<1)
                %% No information about shape factor is provided so we set it to 1
                transferfunction.shape_factor_name = 'ConstantShapeFactor';
                block_dimension = [1,1,1];
                shape_factor_value = 1;
                shape_factor_handle = str2func(transferfunction.shape_factor_name);
                transferfunction.shape_factor_object = shape_factor_handle(block_dimension, shape_factor_value);
            else
                if (nargin<2)
                    msg = ['ERROR: If you provide the shape factor name, make sure to',...
                          'provide also the fracture spacing, [lx,ly,lz]'];
                    error(msg);
                else
                    transferfunction.shape_factor_name = shape_factor_name;
                    shape_factor_handle = str2func(transferfunction.shape_factor_name);
                    transferfunction.shape_factor_object = shape_factor_handle(block_dimension);
                end
            end


        end


        function [Talpha] = calculate_transfer(ktf, model, fracture_fields, matrix_fields)

            %% All calculate_transfer method should have this call. This is a "sanity check" that
            % ensures that the correct structures are being sent to calculate the transfer
            ktf.validate_fracture_matrix_structures(fracture_fields, matrix_fields);
                              

            %% The varibles
            pm = matrix_fields.pm;
            p = fracture_fields.pf;


            %% Additional Properties
            muwm = model.fluid_matrix.muW(pm);
            bWm = model.fluid.bW(pm);

            %% Shape Factor
            %lx,ly,lz are the inverse of the fracture spacing given
            % by the user i.e. the matrix block dimensions
            [sigma]=ktf.shape_factor_object.calculate_shape_factor(model);
            
            %% Compute Transfer
            %(units 1/T)
            kx = model.rock_matrix.perm(:,1);
            try
                ky = model.rock_matrix.perm(:,2);
                kz = model.rock_matrix.perm(:,3);
            catch
                ky = model.rock_matrix.perm(:,1);
                kz = model.rock_matrix.perm(:,1);
            end
            
            if strcmp(ktf.shape_factor_name, 'ConstantShapeFactor')
                sigma = sigma./((1/3)*(kx+ky+kz));
            end
            
            tw=(sigma.*(kx.*ky.*kz).^(1/3).*bWm./muwm).*(p-pm);
            Talpha = tw;
        end

        function [] = validate_fracture_matrix_structures(ktf,fracture_fields, matrix_fields)
            %% We use the superclass to validate the structures of matrix/fracture variables
            validate_fracture_matrix_structures@TransferFunction(ktf,fracture_fields, matrix_fields);
        end

    end


end
