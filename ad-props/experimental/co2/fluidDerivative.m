function u=fluidDerivative(fluid,p,T,property)
% calculate joule tomson coefficent

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

[pa,Ta] = initVariablesADI(double(p),double(T));
switch property
    case 'jouleThomson'
        val=fluid.enthalpy(pa,Ta);
        u=-full(diag(val.jac{1}))./full(diag(val.jac{2}));
    case 'freeExpansion'
        val=fluid.enthalpy(pa,Ta)-pa./fluid.density(pa,Ta);
        u=-full(diag(val.jac{1}))./full(diag(val.jac{2}));
    case 'adiabaticExpansion'
        %val=fluid.enthalpy(pa,Ta)-pa./double(fluid.density(pa,Ta));
        val=fluid.enthalpy(pa,Ta)-pa./fluid.density(pa,Ta)  + double(pa).*(1./fluid.density(pa,Ta));
        u=-full(diag(val.jac{1}))./full(diag(val.jac{2}));
    case 'c_p'
        val=fluid.enthalpy(pa,Ta);
        u=full(val.jac{2});
    case 'comp_p'
        val=fluid.density(pa,Ta);
        u=full(val.jac{1})/double(val);
    case 'comp_t'
        val=fluid.density(pa,Ta);
        u=-full(val.jac{2})./double(val);
    otherwise
        error('no such property implemented')
end
end
