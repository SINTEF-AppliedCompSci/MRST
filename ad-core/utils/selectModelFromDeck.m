function model = selectModelFromDeck(G, rock, fluid, deck, varargin)
    require ad-blackoil

    rs = deck.RUNSPEC;
    check = @(name) isfield(rs, upper(name)) && rs.(upper(name));

    hasgas  = check('gas');
    hasoil  = check('oil');
    haswat  = check('water');
    haspoly = check('polymer');
    
    if haspoly
        assert(~hasgas, 'Polymer model does not support gas');
        model = OilWaterPolymerModel(G, rock, fluid, 'inputdata', deck, varargin{:});
    elseif hasgas && hasoil && haswat
        model = ThreePhaseBlackOilModel(G, rock, fluid, 'inputdata', deck, varargin{:});
    elseif hasoil && haswat
        model = TwoPhaseOilWaterModel(G, rock, fluid, 'inputdata', deck, varargin{:});
    else
        error('Did not find matching model');
    end
end
