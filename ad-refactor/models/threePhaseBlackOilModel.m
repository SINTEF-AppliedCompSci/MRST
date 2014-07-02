classdef threePhaseBlackOilModel < physicalModel
    % Three phase with optional dissolved gas and vaporized oil
    properties
        % Determines if gas can be dissolved into the oil phase
        disgas
        % Determines if oil can be vaporized into the gas phase
        vapoil
        
        % Maximum Rs increment
        drsMax
        
        % CNV style tolerances
        useCNVConvergence
        toleranceCNV;
        toleranceMB;
    end
    
    methods
        function model = threePhaseBlackOilModel(G, rock, fluid, varargin)
            
            model = model@physicalModel(G, rock, fluid);
            
            % Typical black oil is disgas / dead oil, but all combinations
            % are supported
            model.vapoil = false;
            model.disgas = true;
           
            % Max increments
            model.drsMax = inf;
            
            model.useCNVConvergence = true;
            model.toleranceCNV = 1e-7;
            model.toleranceMB = 1e-3;
            
                        
            % All phases are present
            model.oil = true;
            model.gas = true;
            model.water = true;
            
            model = merge_options(model, varargin{:});
            
            d = model.inputdata;
            if ~isempty(d)
                if isfield(d, 'RUNSPEC')
                    if isfield(d, 'VAPOIL')
                        model.vapoil = model.vapoil || d.RUNSPEC.VAPOIL;
                    end
                    if isfield(d.RUNSPEC, 'DISGAS')
                        model.disgas = model.disgas || d.RUNSPEC.DISGAS;
                    end
                else
                    error('Unknown dataset format!')
                end
            end
            model.name = 'BlackOil_3ph';
            model = model.setupOperators(G, rock, 'deck', model.inputdata);
        end
        
        function [convergence, values] = checkConvergence(model, problem, varargin)
            if model.useCNVConvergence
                % Use convergence model similar to commercial simulator
                [convergence, values] = CNV_MBConvergence(model, problem);
            else
                % Use strict tolerances on the residual without any 
                % fingerspitzengefuhlen by calling the parent class
                [convergence, values] = checkConvergence@physicalModel(model, problem, varargin{:});
            end
            
        end
        
        function [problem, state] = getEquations(model, state0, state, dt, drivingForces, varargin)
            [problem, state] = equationsBlackOil(state0, state, dt, ...
                            model.G,...
                            drivingForces,...
                            model.operators,...
                            model.fluid,...
                            'disgas', model.disgas, ...
                            'vapoil', model.vapoil, ...
                            'oil',    model.oil, ...
                            'gas',    model.gas, ...
                            'water',  model.water, ...
                            varargin{:});
            
        end
        
        function [state, report] = updateState(model, state, problem, dx, drivingForces)
            state = updateStateBlackOil(state, dx, drivingForces, model.fluid, ...
             'dsMax',       model.dsMax, ...
             'dpMax',       model.dpMax, ...
             'drsMax',      model.drsMax, ...
             'disgas',      model.disgas, ...
             'vapoil',      model.vapoil);
            report = struct();
        end
    end
end