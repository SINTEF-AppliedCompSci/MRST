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
            opt = struct('deck',                [], ...
                         'drsMax',              inf, ...
                         'dpMax',               inf, ...
                         'dsMax',               .2, ...
                         'disgas',              false,...
                         'useCNVConvergence',   true, ...
                         'toleranceMB',         1e-7, ...
                         'toleranceCNV' ,       1e-3, ...
                         'vapoil',              false);
            opt = merge_options(opt, varargin{:});
            
            model.vapoil = opt.vapoil;
            model.disgas = opt.disgas;
            model.fluid  = fluid;
            model.G   = G;
            
            %
            model.oil = true;
            model.gas = true;
            model.water = true;
            
            % Max increments
            model.drsMax = opt.drsMax;
            model.dpMax = opt.dpMax;
            model.dsMax = opt.dsMax;
            
            model.useCNVConvergence = opt.useCNVConvergence;
            model.toleranceCNV = opt.toleranceCNV;
            model.toleranceMB = opt.toleranceMB;
            
            if ~isempty(opt.deck)
                if isfield(opt.deck.RUNSPEC, 'VAPOIL')
                    model.vapoil = model.vapoil || opt.deck.RUNSPEC.VAPOIL;
                end
                if isfield(opt.deck.RUNSPEC, 'DISGAS')
                    model.disgas = model.disgas || opt.deck.RUNSPEC.DISGAS;
                end
            end
            model.name = 'BlackOil_3ph';
            model = model.setupOperators(G, rock, 'deck', opt.deck);
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