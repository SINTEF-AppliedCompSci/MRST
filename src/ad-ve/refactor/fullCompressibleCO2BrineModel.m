classdef fullCompressibleCO2BrineModel < physicalModel

    % ----------------------------- Class properties -----------------------------
    properties
        
        % Reference temperature at a given reference depth; temperature gradient
        T_ref       % in Kelvin                 
        T_ref_depth % in meters below sea level 
        T_grad      % in degrees Kelvin / 1000 m
        
    end
    % ------------------------------- Class methods -------------------------------
    methods

        %% Constructor
        function model = fullCompressibleCO2BrineModel(Gt, tinfo, varargin)
            opt = struct; % determine which defaults to add here
            opt = merge_options(opt, varargin{:});
            
            % Inherited properties
            model.fluid = ;
            model.G     = Gt;
            model.name  = 'fully_compressible_CO2_brine';
            model       = model.setupOperators(Gt, rock, varargin{:});
            model.oil   = false;
            model.gas   = true;
            model.water = true;
            
            % Other properties
            model.T_ref       = tinfo{1}; 
            model.T_ref_depth = tinfo{2}; 
            model.T_grad      = tinfo{3}; 
          
        end
        
        %% (overridden) member methods
        
        function model = setupOperators(model, G, rock, varargin)
            
        % Since we are working on a top surface grid, we cannot use
        % setupSimComp directly. We modify the permeability of the rock
        % first, in order to make it not the vertical average, but the
        % vertical sum.
            rock.perm = rock.perm .* Gt.cells.H;
            model.operators = setupSimComp(Gt, rock, varargin{:});
        end
        
        function [problem, state] = ...
                getEquations(model, state0, state, dt, drivingForces, varargin)

            tinfo = {model.T_ref, model.T_ref_depth, model.T_grad};
            
            [problem, state] = ...
                equationsCO2BrineCompressible(state0, state, dt, model.G, ...
                                              drivingForces, model.operators, ...
                                              model.fluid, ...
                                              tinfo, ...
                                              varargin{:});
        end
        
        function updateState(model, state, problem, dx, drivingForces)
            
        end
    end
end
