function p = initializeEquilibriumPressures(model, region)
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

    actPh = model.getActivePhases();
    nph = sum(actPh);
    
    cells = region.cells;
    if ischar(cells)
        cells = 1:model.G.cells.num;
    end
    nc = numel(cells);
    
    z = model.G.cells.centroids(cells, 3);

    norm_g = norm(model.gravity);
    if norm_g == 0
        warning('Gravity is zero. Initialization results may not be what you expect.');
    end
    zmax = max(z);
    zmin = min(z);

    zmax = max(zmax, max(region.contacts));
    zmin = min(zmin, min(region.contacts));
    
    p = zeros(nc, nph);
    % Deal with reference phase first
    ix = region.reference_index;

    [p(:, ix), po_fn] = solvePressure(region.datum_pressure, region.datum_depth, region.rho{ix}, norm_g, z, zmin, zmax);

    % Then deal with other phases
    contactNo = 1;
    for ix = 1:nph
        if ix == region.reference_index
            continue
        end
        contact = region.contacts(contactNo);
        
        p_ref = po_fn(contact) + region.pc_sign(ix)*region.contacts_pc(contactNo);
        
        p(:, ix) = solvePressure(p_ref, contact, region.rho{ix}, norm_g, z, zmin, zmax);
        contactNo = contactNo + 1;
    end
end

function [sol_u, sol_d] = integratePressureDrop(dp, P0, Z0, Z_min, Z_max)
    odeopts = odeset('AbsTol', 1.0e-10, 'RelTol', 5.0e-8);
    sol_u = [];
    Z_dist = [Z0, Z_min];
    if ~(Z0 < Z_min) && abs(diff(Z_dist)) > 0
        sol_u = ode45(dp, Z_dist, P0, odeopts);
    end

    sol_d = [];
    Z_dist = [Z0, Z_max];
    if ~(Z0 > Z_max) && abs(diff(Z_dist)) > 0
        sol_d = ode45(dp, Z_dist, P0, odeopts);
    end
end

function [p, p_fn] = solvePressure(P0, Z0, rho, g, Z, Z_min, Z_max)
    dp = @(z, p) g*rho(p, z);
    [sol_u, sol_d] = integratePressureDrop(dp, P0, Z0, Z_min, Z_max);
    
    p = evaluate_pressure(Z, Z0, P0, sol_u, sol_d);
    
    p_fn = @(z) evaluate_pressure(z, Z0, P0, sol_u, sol_d);
end

function p = evaluate_pressure(Z, Z0, P0, sol_u, sol_d)
	p = zeros(size(Z));
    down = Z > Z0;
    if any(down)
        p(down) = deval(sol_d, Z(down));
    end
    
    above = Z < Z0;
    if any(above)
        p(above) = deval(sol_u, Z(above));
    end
    p(Z == Z0) = P0;
end
