function initState = computeInitAquifer(model, initState, initaqvolumes)
%Undocumented Utility Function

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

    aquifers = model.aquifers;
    aquind   = model.aquind;
    aquid    = aquifers(:, aquind.aquid);
    naq      = max(aquid);
    nconn    = size(aquid, 1);
    
    p_aq   = zeros(naq, 1);
    p_aqAD = initVariablesADI(p_aq);
    
    state = initState;
    state = model.setProp(state, 'aquiferpressures', p_aqAD);
    state = model.setProp(state, 'aquifervolumes'  , initaqvolumes);
    
    dt = 0;
    
    q = computeAquiferFluxes(model, state, dt);
    % We want to compute p_aq that minimizes q'*q. The function q is affine with respect to p_aq,
    % that is, q = M*p_aq + r, for some M and r (see function
    % computeAquiferFluxes). The solution p_aq of the minimization problem is
    % then p_aq = -inv(M'*M)*(M'*r)
    r = q.val;
    M = q.jac{1};
    
    p_aq = -(M'*M)\(M'*r);
    initState = model.setProp(initState, 'aquiferpressures', p_aq);
    initState = model.setProp(initState, 'aquifervolumes', V_aq);
end
