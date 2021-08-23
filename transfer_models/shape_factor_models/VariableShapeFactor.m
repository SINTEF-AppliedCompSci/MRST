classdef VariableShapeFactor < ShapeFactor
    % Class representing a shape factor which is potentially different 
    % for every cell.    
    % 
    % INITIALIZATION:
    %    Initialized from a transfer function constructor
    %
    % METHODS:
    %   calculate_shape_factor - Return the shape factors for all cells
    
    % ============================= Class properties ==========================
    properties
        smin   % The (linear) shape factor
    end

    % ============================== Public methods ===========================
    methods

        % Constructor
        function shapefactor = VariableShapeFactor(par) %, qmax, pow)
            
            % Interpret the input parameters as [rock.smin block_dimension]
            smin = par(:, 1);
            block_dimension = par(:, 2:4);
            
            % Call superclass constructor to set the block_dimension
            shapefactor = shapefactor@ShapeFactor(block_dimension);

            % Initialize class properties 
            shapefactor.smin = smin;
            
        end
        
        % Return the shape factors times matrix permeability for all cells
        function [sigma] = calculate_shape_factor(sf, model)      
                    
            % Use the average matrix permeability
            kx = model.rock_matrix.perm(:,1);
            try
                ky = model.rock_matrix.perm(:,2);
                kz = model.rock_matrix.perm(:,3);
            catch
                ky = model.rock_matrix.perm(:,1);
                kz = model.rock_matrix.perm(:,1);
            end      
            k = (kx + ky + kz) / 3;
            
            sigma = sf.smin .* k;
          
        end
        
    end
end
