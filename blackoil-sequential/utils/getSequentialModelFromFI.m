function model = getSequentialModelFromFI(fimodel, varargin)
% For a given fully implicit model, output the corresponding pressure/transport model 
    rock  = fimodel.rock;
    fluid = fimodel.fluid;
    G     = fimodel.G;
    
    switch lower(class(fimodel))
        case 'twophaseoilwatermodel'
            pressureModel  = PressureOilWaterModel(G, rock, fluid);
            transportModel = TransportOilWaterModel(G, rock, fluid);
        case 'threephaseblackoilmodel'
            pressureModel  = PressureBlackOilModel(G, rock, fluid);
            transportModel = TransportBlackOilModel(G, rock, fluid);
        otherwise
            error('mrst:getSequentialModelFromFI', ...
            ['Sequential model not implemented for ''' class(fimodel), '''']);
    end
    model = SequentialPressureTransportModel(pressureModel, transportModel, varargin{:});
end
