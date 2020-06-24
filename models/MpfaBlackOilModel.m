classdef MpfaBlackOilModel < GenericBlackOilModel

    properties

        eta
        bcetazero
        
    end

    methods
        
        function model = MpfaBlackOilModel(G, rock, fluid, varargin)
            
            model = model@GenericBlackOilModel(G, rock, fluid, varargin{:});
            model = merge_options(model, varargin{:});
            
            % Process the grid for mechanical computation
            if ~ismember('createAugmentedGrid', model.G.type)
                model.G = createAugmentedGrid(model.G);
            end
        
            model.eta = 0;
            model.bcetazero = false;
            
            % Add mechanical operators  
            model.operators = setupMpfaAdOperators(model);

        end


        function model = setupStateFunctionGroupings(model, varargin) 
            
            model = setupStateFunctionGroupings@GenericBlackOilModel(model, varargin{:});
            fd = model.FluxDiscretization;
            fd = fd.setStateFunction('PermeabilityPotentialGradient', MpfaKgrad(model));
            model.FluxDiscretization = fd;
            
        end


    end

end
