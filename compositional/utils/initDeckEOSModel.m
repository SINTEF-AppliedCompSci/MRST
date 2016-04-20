function model = initDeckEOSModel(deck)
% Set up a EOS model from a parsed deck
    names = deck.PROPS.CNAMES;
    Tcrit = deck.PROPS.TCRIT';
    Pcrit = deck.PROPS.PCRIT';
    Vcrit = deck.PROPS.VCRIT';
    acf   = deck.PROPS.ACF';
    mw    = deck.PROPS.MW';
    
    T_ref = 273.15 + 15;
    P_ref = 1*atm;
    
    fluid = CompositionalFluid(names, Tcrit, Pcrit, Vcrit, acf, mw, T_ref, P_ref);
    model = EquationOfStateModel([], fluid);
end