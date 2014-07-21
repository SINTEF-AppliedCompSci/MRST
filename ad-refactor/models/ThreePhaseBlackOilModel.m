classdef ThreePhaseBlackOilModel < PhysicalModel
    % Three phase with optional dissolved gas and vaporized oil
    properties
        % Determines if gas can be dissolved into the oil phase
        disgas
        % Determines if oil can be vaporized into the gas phase
        vapoil
        
        % Maximum Rs/Rv increment
        drsMax
        
        % Use alternate tolerance scheme
        useCNVConvergence
        
        % CNV tolerance (inf-norm-like)
        toleranceCNV;
        
        % MB tolerance values (2-norm-like)
        toleranceMB;
        % Well tolerance if CNV is being used
        toleranceWellBHP;
        % Well tolerance if CNV is being used
        toleranceWellRate;
        
        % Update wells
        
    end
    
    methods
        function model = ThreePhaseBlackOilModel(G, rock, fluid, varargin)
            
            model = model@PhysicalModel(G, rock, fluid);
            
            % Typical black oil is disgas / dead oil, but all combinations
            % are supported
            model.vapoil = false;
            model.disgas = true;
           
            % Max increments
            model.drsMax = inf;
            
            model.useCNVConvergence = true;
            model.toleranceCNV = 1e-3;
            model.toleranceMB = 1e-7;
            model.toleranceWellBHP = 1*barsa;
            model.toleranceWellRate = 1/day;
                        
            % All phases are present
            model.oil = true;
            model.gas = true;
            model.water = true;
            model.componentNames = {'sw', 'so', 'sg'};
            
            model = merge_options(model, varargin{:});
            
            d = model.inputdata;
            if ~isempty(d)
                if isfield(d, 'RUNSPEC')
                    if isfield(d.RUNSPEC, 'VAPOIL')
                        model.vapoil = d.RUNSPEC.VAPOIL;
                    end
                    if isfield(d.RUNSPEC, 'DISGAS')
                        model.disgas = d.RUNSPEC.DISGAS;
                    end
                else
                    error('Unknown dataset format!')
                end
            end
            model.name = 'BlackOil_3ph';
            model = model.setupOperators(G, rock, 'deck', model.inputdata);
        end
        
        function [fn, index] = getVariableField(model, name)
            switch(lower(name))
                case 'rs'
                    fn = 'rs';
                    index = 1;
                case 'rv'
                    fn = 'rv';
                    index = 1;
                otherwise
                    % Basic phases are known to the base class
                    [fn, index] = getVariableField@PhysicalModel(model, name);
            end
        end

        function [convergence, values] = checkConvergence(model, problem, varargin)
            if model.useCNVConvergence
                % Use convergence model similar to commercial simulator
                [conv_cells, v_cells] = CNV_MBConvergence(model, problem);
                [conv_wells, v_wells] = checkWellConvergence(model, problem);
                
                convergence = all(conv_cells) && all(conv_wells);
                values = [v_cells, v_wells];
            else
                % Use strict tolerances on the residual without any 
                % fingerspitzengefuhlen by calling the parent class
                [convergence, values] = checkConvergence@PhysicalModel(model, problem, varargin{:});
            end            
        end
        
        function [problem, state] = getEquations(model, state0, state, dt, drivingForces, varargin)
            [problem, state] = equationsBlackOil(state0, state, model, dt, ...
                            drivingForces, varargin{:});
            
        end
        
        function [state, report] = updateState(model, state, problem, dx, drivingForces)
            state = updateStateBlackOilGeneric(model, state, problem, dx, drivingForces);
            report = struct();
        end
    end
end