classdef GenericTracerModelDG < TransportModelDG
    
    properties
    end
    methods
        function model = GenericTracerModelDG(pmodel, varargin)
            fluid = initSimpleADIFluid('phases', 'WO', 'n', [1,1], 'rho', [1,1], 'mu', [1,1]);
            pmodel.fluid = fluid;
            pmodel.water = true;
            pmodel.oil   = true;
            pmodel.gas   = false;
            pmodel.useCNVConvergence = false;
            pmodel.nonlinearTolerance = 1e-3;
            pmodel.Components{1} = ImmiscibleComponent('Tracer', 1);
            pmodel.Components{2} = ImmiscibleComponent('Dummy' , 2);
            model = model@TransportModelDG(pmodel, 'formulation', 'missingPhase', varargin{:});
        end
    end
    
end