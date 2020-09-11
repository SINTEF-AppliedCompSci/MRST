classdef ConstantShapeFactor < ShapeFactor
    %CONSTANT SHAPE FACTOR
    %
    % class to provide a constant given shape factor

    properties
        shape_factor_value;
    end

    methods
        function shapefactor = ConstantShapeFactor(block_dimension, shape_factor_value)
            has_permeability = 1;
            shapefactor = shapefactor@ShapeFactor(block_dimension, has_permeability);
            if nargin<4
              shapefactor.shape_factor_value = 1;
            else
              shapefactor.shape_factor_value = shape_factor_value;
            end
        end

        function [sigma] = calculateShapeFactor(sf, rock)
            
            kx = rock.perm(:,1);
            try
                ky = rock.perm(:,2);
                kz = rock.perm(:,3);
            catch
                ky = rock.perm(:,1);
                kz = rock.perm(:,1);
            end

            kavg = (1/3)*(kx+ky+kz);
            sigma = kavg.*sf.shape_factor_value;

        end
    end
end
