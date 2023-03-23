classdef KazemiSinglePhaseTransferFunction < TransferFunction
 
    properties
        shape_factor_object
        shape_factor_name
    end
    
    methods
        function transferfunction = KazemiSinglePhaseTransferFunction(shape_factor_name,fracture_spacing)
            
            transferfunction = transferfunction@TransferFunction();
            transferfunction.nphases = 1;
            
            if (nargin<1)
                %% No information about shape factor is provided so we set it to 1
                transferfunction.shape_factor_name = 'ConstantShapeFactor';
                fracture_spacing = [1,1,1];
                shape_factor_value = 1;
                shape_factor_handle = str2func(transferfunction.shape_factor_name);
                transferfunction.shape_factor_object = shape_factor_handle(fracture_spacing,shape_factor_value);
            else
                if (nargin<2)
                    msg = ['ERROR: If you provide the shape factor name, make sure to',...
                          'provide also the fracture spacing, [lx,ly,lz]'];
                    error(msg);
                else
                    transferfunction.shape_factor_name = shape_factor_name;
                    shape_factor_handle = str2func(transferfunction.shape_factor_name);
                    transferfunction.shape_factor_object = shape_factor_handle(fracture_spacing);
                end
            end
            
        end
        
        function [Talpha] = calculate_transfer(ktf,model,fracture_fields,matrix_fields)

            %% All calculate_transfer method should have this call. This is a "sanity check" that
            % ensures that the correct structures are being sent to calculate the transfer
            ktf.validate_fracture_matrix_structures(fracture_fields,matrix_fields);                                         
                                                                                  
            %% The varibles
            pm = matrix_fields.pm;
            pf = fracture_fields.pf;
            
            %% Additional Properties
            muwm = model.fluid_matrix.muW(pm);
            bWm = model.fluid_matrix.bW(pm);
			
			%% Shape Factor
            [sigma]=ktf.shape_factor_object.calculate_shape_factor(model);
           
            %% Compute Transfer
            %(units 1/T)
            % Note: shape factors include permeability
            tw=(sigma.*bWm./muwm).*(pf-pm); 
            
            %% Note that we return a 2x1 Transfer since our model is 2ph
            Talpha{1} = tw;
        end
        
        function [] = validate_fracture_matrix_structures(ktf,fracture_fields,matrix_fields)
            %% We use the superclass to validate the structures of matrix/fracture variables                                          
            validate_fracture_matrix_structures@TransferFunction(ktf,fracture_fields,matrix_fields);
        end
        
    end
    
    
end
