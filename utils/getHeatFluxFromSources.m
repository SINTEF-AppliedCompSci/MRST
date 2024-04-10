function forces = getHeatFluxFromSources(model, forces, drivingForces)
%Add heat flux from boundary conditions

    % Early return if no BCs are given
    [forces.bc.heatFlux, forces.src.heatFlux] = deal([]);
    if ~isempty(forces.bc.sourceCells)
        bc = drivingForces.bc;
        propsRes = bc.propsRes;
        propsBC  = bc.propsForce;
        [qAdv, q] = computeAdvectiveHeatFlux( ...
            model, propsRes, propsBC, forces.bc.phaseMass);
        qCond = computeConductiveHeatFlux(propsRes, propsBC, bc);
        forces.bc.advHeatFlux = q;
        forces.bc.condHeatFlux = qCond;
        forces.bc.heatFlux = qAdv + qCond;
    end
    
    if ~isempty(forces.src.sourceCells)
        src = drivingForces.src;
        propsRes  = src.propsRes;
        propsSrc  = src.propsForce;
        [qAdv, q] = computeAdvectiveHeatFlux( ...
            model, propsRes, propsSrc, forces.src.phaseMass);
        qCond = computeConductiveHeatFlux(propsRes, propsSrc, src);
        forces.src.advHeatFlux = q;
        forces.src.condHeatFlux = qCond;
        forces.src.heatFlux = qAdv + qCond;
    end
    
end

%-------------------------------------------------------------------------%
function [q, qph] = computeAdvectiveHeatFlux(model, propsRes, propsBC, src)

    nph = model.getNumberOfPhases();
    q   = 0;
    qph = cell(1,nph);
    hbc = model.getProps(propsBC , 'PhaseEnthalpy');
    hr  = model.getProps(propsRes, 'PhaseEnthalpy');    
    for i = 1:nph
        inflow = src{i} > 0;
        h      = inflow.*hbc{i} + ~inflow.*hr{i};
        qph{i} = src{i}.*h;
        q      = q + qph{i};
    end
    
end

%-------------------------------------------------------------------------%
function q = computeConductiveHeatFlux(propsRes, propsBC, bc)

    is_Hflux = ~isnan(bc.Hflux);
    q = -(propsRes.Thr + propsRes.Thf).*(propsRes.T - propsBC.T);
    q(is_Hflux) = bc.Hflux(is_Hflux);
    
end

%{
Copyright 2009-2024 SINTEF Digital, Mathematics & Cybernetics.

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