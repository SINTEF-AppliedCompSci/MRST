function operators = setupMpfaOperators(model)
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
    bcstruct = model.fluid.bcstruct;

    eta = 0;
    bcetazero = false;

    perm = model.rock.perm;

    assert(numel(perm) == G.cells.num, 'only isotropic perm for the moment');
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
    src = []; % no source at this stage

    assembly = assembleMPFA(G, K, bcstruct, src, eta, tbls, mappings, 'addAdOperators', true);

    operators = assembly.adoperators;

    % recover pore volume from TPFA (could be optimized!)
    rock = model.rock;
    tpfaoperators = setupOperatorsTPFA(G, rock);
    operators.pv = tpfaoperators.pv;
end
