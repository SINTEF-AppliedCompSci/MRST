classdef ShapeFactor
    %SHAPE_FACTOR Summary of this class goes here
    %   Detailed explanation goes here
    %PROPERTIES
    %lx,ly,lz are the dimensions of the matrix block
    % these are claculated as the inverse of the fracture spacing
    % provided by the user

    properties
        block_dimension
    end

    methods
        function shapefactor = ShapeFactor(block_dimension)
            shapefactor.block_dimension = block_dimension;
        end

        function [sigma] = calculate_shape_factor(sf,model)
            sigma=1;
        end
    end
end
