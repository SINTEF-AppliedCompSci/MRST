function model = initDeckEOSModel(deck)
% Set up a EOS model from a parsed deck
    names = deck.PROPS.CNAMES;
    Tcrit = deck.PROPS.TCRIT';
    Pcrit = deck.PROPS.PCRIT';
    if isfield(deck.PROPS, 'VCRIT')
        Vcrit = deck.PROPS.VCRIT';
    else
        R = 8.3144598;
        Vcrit = R*deck.PROPS.ZCRIT'.*Tcrit./Pcrit;
    end
    acf   = deck.PROPS.ACF';
    mw    = deck.PROPS.MW';
    
    T_ref = 273.15 + 15;
    P_ref = 1*atm;
    
    fluid = CompositionalFluid(names, Tcrit, Pcrit, Vcrit, acf, mw, T_ref, P_ref);
    model = EquationOfStateModel([], fluid);
end

%{
Copyright 2009-2017 SINTEF ICT, Applied Mathematics.

This file is part of The MATLAB Reservoir Simulation Toolbox (MRST).

MRST is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

MRST is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with MRST.  If not, see <http://www.gnu.org/licenses/>.
%}
