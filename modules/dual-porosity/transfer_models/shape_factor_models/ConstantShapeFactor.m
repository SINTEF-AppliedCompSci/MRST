classdef ConstantShapeFactor < ShapeFactor
    %CONSTANT SHAPE FACTOR
    %
    % class to provide a constant given shape factor

    properties
        shape_factor_value;
    end

    methods
        function shapefactor = ConstantShapeFactor(block_dimension,shape_factor_value)
          shapefactor = shapefactor@ShapeFactor(block_dimension);
          if nargin<4
              shapefactor.shape_factor_value = 1;
          else
              shapefactor.shape_factor_value = shape_factor_value;
          end
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

            kavg = (1/3)*(kx+ky+kz);
            sigma = kavg.*sf.shape_factor_value;

        end
    end
end
