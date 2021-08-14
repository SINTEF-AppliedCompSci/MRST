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
            
            % Call superclass constructor to satisfy the syntax of EclipseTransferFunction
            shapefactor = shapefactor@ShapeFactor(block_dimension);

            % Initialize class properties 
            shapefactor.smin = smin;
            
        end
        
%                 % Constructor
%         function shapefactor = VariableShapeFactor(smin) %, qmax, pow)
%             
%             % Initialize class properties 
%             shapefactor.smin = smin;
%             
%         end
        
        
        

        % Return the shape factors times matrix permeability for all cells
        function [sigma] = calculate_shape_factor(sf, model)      
                    
            % Use the average matrix permeability
            kxx = model.rock_matrix.perm(:,1);
            try
                kyy = model.rock_matrix.perm(:,3);
            catch
                kyy = model.rock_matrix.perm(:,1);
            end
            k = 0.5*(kxx + kyy);
            
            sigma = sf.smin .* k;
        end
        
    end
end
