function frac = getNonUnitMassFraction(model, molfraction, molfraction_normalized)
% Internal utility. Intentionally undocumented.

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

    model = model.EOSModel;
    ncomp = numel(molfraction);
    mass = cell(1, ncomp);
    total = 0;
    mw = model.fluid.molarMass;
    if iscell(molfraction)
        for i = 1:numel(molfraction)
            mass{i} = molfraction{i}.*mw(i);
            total = total + molfraction_normalized{i}.*mw(i);
        end
        frac = cellfun(@(x) x./total, mass, 'UniformOutput', false);
    else
        mass = molfraction.*mw;
        total = sum(molfraction_normalized.*mw, 2);
        frac = mass./total;
    end
    
end
