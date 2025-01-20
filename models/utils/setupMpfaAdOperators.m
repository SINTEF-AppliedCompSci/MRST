function operators = setupMpfaAdOperators(model, varargin)    
% We setup the usual operators and the mpfa flux operator
%
% we set up the mappings
%
%   F1 : nodefacetbl -> intfacetbl
%   F2 : celltbl -> intfacetbl
%
%  such that the flux u in intfacetbl (interior faces) is given by
%
%  u  = [F1  F2] * [ pnf (pressure at nodefacetbl);
%                    pc  (pressure at celltbl)];
%
%  Then, we proceed with the reduction to remove dependency in pnf
%  (pressure at nodefacetbl)

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


    opt = struct('useVirtual', false);
    opt = merge_options(opt, varargin{:});

    useVirtual = opt.useVirtual;
    
    G         = model.G;
    rock      = model.rock;
    eta       = model.eta;
    bcetazero = model.bcetazero;

    perm = rock.perm;

    assert(numel(perm) == G.cells.num, 'only isotropic perm for the moment');

    [tbls, mappings] = setupMpfaStandardTables(G, 'useVirtual', false);

    celltbl      = tbls.celltbl;
    vec12tbl     = tbls.vec12tbl;
    cellvec12tbl = tbls.cellvec12tbl;

    prod = TensorProd();
    prod.tbl1 = vec12tbl;
    prod.tbl2 = celltbl;
    prod.tbl3 = cellvec12tbl;

    if useVirtual
        
        prod.pivottbl = cellvec12tbl;

        N = tbls.celltbl.num;
        [vec2, vec1, i] = ind2sub([d_num, d_num, N], (1 : cellvec12tbl.num)');

        prod.dispind1 = sub2ind([d_num, d_num], vec2, vec1);

        prod.dispind2 = i;

        prod.dispind3 = (1 : cellvec12tbl.num)';
        
        prod.issetup = true;
        
    else
        prod = prod.setup();
    end

    K = prod.eval([1; 0; 0; 1], perm);
    src = []; % no source at this stage
    bcstruct.bcdirichlet = []; % no Dirichlet
    bcstruct.bcneumann = []; % zero neumann bc

    assembly = assembleMPFA(G, K, bcstruct, src, eta, tbls, mappings, 'addAdOperators', true);

    operators = setupOperatorsTPFA(G, rock);
    operators.fluxop = assembly.adoperators.fluxop;
end
