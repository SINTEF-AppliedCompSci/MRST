classdef threePhaseBlackOilModel < physicalModel
    % Three phase with optional dissolved gas and vaporized oil
    properties
        disgas
        vapoil
        
        % Maximum rs increment
        drsMax
    end
    
    methods
        function model = threePhaseBlackOilModel(G, rock, fluid, varargin)
            opt = struct('deck',                [], ...
                         'disgas',              false,...
                         'vapoil',              false);
            opt = merge_options(opt, varargin{:});
            
            model.vapoil = opt.vapoil;
            model.disgas = opt.disgas;
            model.fluid  = fluid;
            model.G   = G;
            
            % Max increments
            model.drsMax = inf;
            
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
        
        function problem = getEquations(model, state0, state, dt, drivingForces, varargin)
            problem = equationsBlackOil(state0, state, dt, ...
                            model.G,...
                            drivingForces,...
                            model.operators,...
                            model.fluid,...
                            'disgas', model.disgas, ...
                            'vapoil', model.vapoil, ...
                            varargin{:});
            
        end
        
        function state = updateState(model, state, dx, drivingForces)
            state = updateStateBlackOil(state, dx, drivingForces, model.fluid, ...
             'dsMax',       model.dsMax, ...
             'dpMax',       model.dpMax, ...
             'drsMax',      model.drsMax, ...
             'disgas',      model.disgas, ...
             'vapoil',      model.vapoil);
        end
    end
end