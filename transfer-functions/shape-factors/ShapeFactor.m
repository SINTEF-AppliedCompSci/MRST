classdef ShapeFactor
    %SHAPE_FACTOR 
    %lx,ly,lz are the dimensions of the matrix block
    % these are calculated as the inverse of the fracture spacing
    % provided by the user

    properties
        block_dimension
        has_permeability
    end

    methods
        function shapefactor = ShapeFactor(block_dimension, has_permeability)
            shapefactor.block_dimension = block_dimension;
            shapefactor.has_permeability = has_permeability;
        end

        function [sigma] = calculateShapeFactor(sf, rock)
            sigma = 1;
        end
    end
end
