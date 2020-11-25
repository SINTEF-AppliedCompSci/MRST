function operators = setupBiotOperators(model)
%Undocumented Utility Function

%{
Copyright 2020 University of Bergen and SINTEF Digital, Mathematics & Cybernetics.

This file is part of the MPSA-W module for the MATLAB Reservoir Simulation Toolbox (MRST).

The MPSA-W module is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

The MPSA-W module is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with the MPSA-W module.  If not, see <http://www.gnu.org/licenses/>.
%}


    G = model.G;
    mech  = model.mech;
    fluid = model.fluid;
    rock  = model.rock;

    % setup pore volume
    pv = poreVolume(G, rock);

    % setup mechanic input
    mechprops = mech.prop;
    loadstruct = mech.loadstruct;

    % setup flow input
    perm = model.rock.perm;

    assert(numel(perm) == G.cells.num, 'only isotropic perm for the moment (for simplicity)');
    [tbls, mappings] = setupStandardTables(G, 'useVirtual', false);

    celltbl = tbls.celltbl;
    colrowtbl = tbls.colrowtbl;
    cellcolrowtbl = tbls.cellcolrowtbl;

    prod = TensorProd();
    prod.tbl1 = colrowtbl;
    prod.tbl2 = celltbl;
    prod.tbl3 = cellcolrowtbl;
    prod = prod.setup();

    K = prod.eval([1; 0; 0; 1], perm);

    if isfield(fluid, 'bcstruct')
        bcstruct = fluid.bcstruct;
    else
        bcstruct.bcdirichlet = [];
        bcstruct.bcneumann = [];
        bcstruct.bcneumann = [];
    end

    eta = model.eta;
    bcetazero = model.bcetazero;

    % run assembly
    fluidprops.K = K;
    fluidforces.bcstruct = bcstruct;
    fluidforces.src = []; % no explicit sources here (for the moment)

    coupprops.rho = 0; % the accumulation term is set within the equation
    coupprops.alpha = rock.alpha;
    props = struct('mechprops' , mechprops , ...
                   'fluidprops', fluidprops, ...
                   'coupprops' , coupprops);

    drivingforces = struct('mechanics', loadstruct, ...
                           'fluid'    , fluidforces);

    [tbls, mappings] = setupStandardTables(G);

    assembly = assembleBiot(G, props, drivingforces, eta, tbls, mappings, 'addAdOperators', true);

    operators = assembly.adoperators;
    operators.pv = pv;
    operators.tbls = tbls;
end
