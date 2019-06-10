classdef TransferFunction
    %TRANSFERFN Summary of this class goes here
    %   Detailed explanation goes here
    % DESCRIPTION:
    %   Base Transfer function class
    %

    properties
        nphases
    end
    methods
        function transferfunction = TransferFunction()
            
        end
        
        function [] = validate_fracture_matrix_structures(tf,fracture_fields,matrix_fields)

            fracture_fields = fieldnames(fracture_fields);
            matrix_fields = fieldnames(matrix_fields);
            
            if(tf.nphases == 1)
                assert(any(ismember(fracture_fields,'pf')),...
                      'Could not find "pof" field in "fracture_fields" structure')
                  
                assert(any(ismember(matrix_fields,'pm')),...
                    'Could not find "pom" field in "matrix_fields" structure')
            end
            if(tf.nphases == 2)
                assert(any(ismember(fracture_fields,'pof')),...
                    'Could not find "pof" field in "fracture_fields" structure')
                assert(any(ismember(fracture_fields,'swf')),...
                    'Could not find "swf" field in "fracture_fields" structure')
                
                assert(any(ismember(matrix_fields,'pom')),...
                    'Could not find "pom" field in "matrix_fields" structure')
                assert(any(ismember(matrix_fields,'swm')),...
                    'Could not find "swm" field in "matrix_fields" structure')
            end
            if(tf.nphases == 3) % 3 phases
                assert(any(ismember(fracture_fields,'pof')),...
                    'Could not find "pof" field in "fracture_fields" structure')
                assert(any(ismember(fracture_fields,'swf')),...
                    'Could not find "swf" field in "fracture_fields" structure')
                assert(any(ismember(fracture_fields,'sgf')),...
                    'Could not find "sgf" field in "fracture_fields" structure')
                
                assert(any(ismember(matrix_fields,'pom')),...
                    'Could not find "pom" field in "matrix_fields" structure')
                assert(any(ismember(matrix_fields,'swm')),...
                    'Could not find "swm" field in "matrix_fields" structure')
                assert(any(ismember(matrix_fields,'sgm')),...
                    'Could not find "sgf" field in "matrix_fields" structure')
            end
            
        end
        
        function [Talpha] = calculate_transfer(tf,model,fracture_fields,matrix_fields)
            % This method should be reimplemented in the derived classes
            Talpha{1} = 0;
            Talpha{2} = 0;
            Talpha{3} = 0;
        end
        
    end
    
    
end


%{
Copyright 2009-2016 SINTEF ICT, Applied Mathematics.

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
