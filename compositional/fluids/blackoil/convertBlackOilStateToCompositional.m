function state0 = convertBlackOilStateToCompositional(bomodel, state)
% Convert BO state to compositional-like state

%{
Copyright 2009-2021 SINTEF Digital, Mathematics & Cybernetics.

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

    [p, sO, sG, sW] = bomodel.getProps(state, 'pressure', 'sO', 'sG', 'sW');
    if bomodel.disgas
        rs = bomodel.getProp(state, 'rs');
    else
        rs = 0;
    end
    if bomodel.vapoil
        rv = bomodel.getProp(state, 'rv');
    else
        rv = 0;
    end
    [xo, xg, yo, yg, zo, zg, rhoO, rhoG] = blackOilToMassFraction(bomodel, p, sO, sG, rs, rv);
    
    state0 = initResSol(bomodel.G, p, [sW, sO, sG]);
    state0.components = [zo, zg];
    state0.x = [xo, xg];
    state0.L = rhoO.*sO./(rhoO.*sO + rhoG.*sG);
    
    state0.y = repmat([yo, yg], bomodel.G.cells.num, 1);
    state0.T = repmat(273.15, bomodel.G.cells.num, 1);
end
