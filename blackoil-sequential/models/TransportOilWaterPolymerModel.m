classdef TransportOilWaterPolymerModel < OilWaterPolymerModel
    % Two phase oil/water system with polymer
    properties
        conserveWater
        conserveOil
        staticUpwind
    end
    
    methods
        function model = TransportOilWaterPolymerModel(G, rock, fluid, varargin)
            
            model = model@OilWaterPolymerModel(G, rock, fluid);
            
            model.conserveWater = true;
            model.conserveOil   = false;
            
            model.staticUpwind  = false;

            model = merge_options(model, varargin{:});
            
            assert(~(model.conserveWater && model.conserveOil), ... 
                            'Sequential form only conserves n-1 phases');
            
            % Ensure simple tolerances
            model.useCNVConvergence = false;
        end
        
        function [problem, state] = getEquations(model, state0, state, ...
                dt, drivingForces, varargin)
            [problem, state] = transportEquationOilWaterPolymer(...
                state0, state, model, dt, drivingForces,...
                'solveForOil',   model.conserveOil, ...
                'solveForWater', model.conserveWater, ...
                varargin{:});
        end
        
        function [convergence, values] = checkConvergence(model, problem, varargin)
%             [convergence, values] = checkConvergence@PhysicalModel(model, problem, varargin{:});
            
%             model.oil = false;
%             [convergence, values] = checkConvergence@ReservoirModel(...
%                 model, problem, varargin{:});
%             model.oil = true;

            if model.useCNVConvergence
                
                % Use convergence model similar to commercial simulator
                model.oil = false;
                [conv_cells, v_cells] = CNV_MBConvergence(model, problem);
                model.oil = true;
                
                % We do not include the polymer equation in the convergence
                convergence = all(conv_cells);
                
                values = v_cells;
            else
                % Use strict tolerances on the residual without any 
                % fingerspitzengefuhlen by calling the parent class
                [convergence, values] = checkConvergence@PhysicalModel(...
                    model, problem, varargin{:});
            end  


            % Remove polymer residual from conv criterion
%             values = values(~strcmpi(problem.equationNames, 'polymer'));
%             convergence = all(values < model.nonlinearTolerance);
            
            % Always make at least one update so that the problem actually changes.
            convergence = convergence && problem.iterationNo > 1;
        end
    end
end
