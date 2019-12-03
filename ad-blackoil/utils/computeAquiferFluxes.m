function q = computeAquiferFluxes(model, state, dt)
%Undocumented Utility Function

%{
Copyright 2009-2019 SINTEF Digital, Mathematics & Cybernetics.

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
    
    [p, sW] = model.getProps(state, 'pressure', 'water');
    p_aq = model.getProp(state, 'aquiferpressures');
    V_aq = model.getProp(state, 'aquifervolumes');

    aquifers  = model.aquifers;
    aquind  = model.aquind;
    
    alpha     = aquifers(:, aquind.alpha);
    J         = aquifers(:, aquind.J);
    conn      = aquifers(:, aquind.conn);
    depthconn = aquifers(:, aquind.depthconn);
    depthaq   = aquifers(:, aquind.depthaq);
    C         = aquifers(:, aquind.C);
    aquid     = aquifers(:, aquind.aquid);
    
    nconn = size(conn, 1);
    naq = max(aquid);
    aquid2conn = sparse(aquid, (1 : nconn)', 1, naq, nconn)';
    
    p_aq = aquid2conn*p_aq;
    V_aq = aquid2conn*V_aq;

    p  = p(conn);
    sW = sW(conn);
    
    fluid = model.fluid;
    pcOW = 0;
    ph = model.getPhaseNames();
    isw = (ph == 'W');
    if isfield(fluid, 'pcOW') && ~isempty(sW)
        pcOW = model.getProps(state, 'CapillaryPressure');
        pcOW = pcOW{isw}(conn);
    end

    b = model.getProps(state, 'ShrinkageFactors');
    bW = b{isw}(conn);
    rhoW = bW.*fluid.rhoWS;
    
    Tc = C.*V_aq./J;
    if dt == 0
        coef = 1;
    else
        coef = (1 - exp(-dt./Tc))./(dt./Tc);
    end
    
    g = model.gravity(3);
    q = alpha.*J.*(p_aq + pcOW - p + rhoW.*g.*(depthconn - depthaq)).*coef;
    
end
