function assembly = assembleMPFA(G, K, bcstruct, src, eta, tbls, mappings, varargin)
%
%   
% Before nodal reduction, the solution is given by the system
%
% A = [[A11, A12, -D];
%      [A21, A22,  0];
%      [D' , 0  ,  0]];
%
% u = [pnf     (pressure at nodefacetbl);
%      pc      (pressure at celltbl);
%      lagmult (flux values Dirichlet boundary)];
%
% f = [flux   (flux at nodefacetbl);
%      src    (source at celltbl);
%      bcvals (pressure values at Dirichlet boundary)];
%
% A*u = f
%
% Note: extforce is sparse and should only give contribution at facets
% that are at the boundary
%
% By construction of the method, the matrix A11 is block-diagonal. Hence,
% we invert it directly and reduce to a cell-centered scheme.
% 
%
% First, we assemble the matrices A11, A12, A21, A22

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


    opt = struct('bcetazero'           , true , ...
                 'addAssemblyMatrices' , false, ...
                 'addAdOperators'      , false, ...
                 'onlyAssemblyMatrices', false, ...
                 'useVirtual'          , false);
    opt = merge_options(opt, varargin{:});
    
    useVirtual = opt.useVirtual;
    
    cellvec12tbl       = tbls.cellvec12tbl;
    cellnodevec12tbl   = tbls.cellnodevec12tbl;
    cellnodefacetbl    = tbls.cellnodefacetbl;
    cellnodefacevectbl = tbls.cellnodefacevectbl;
    celltbl            = tbls.celltbl;
    vectbl             = tbls.vectbl;
    nodefacevectbl     = tbls.nodefacevectbl;
    nodefacetbl        = tbls.nodefacetbl;

    % Some shortcuts
    c_num     = celltbl.num;
    cnf_num   = cellnodefacetbl.num;
    nf_num    = nodefacetbl.num;
    d_num     = vectbl.num;
    
    % Setup main assembly matrices
    bcdirichlet = bcstruct.bcdirichlet;
    opts = struct('eta', eta, ...
                  'bcetazero', opt.bcetazero);
    output = coreMpfaAssembly(G, K, bcdirichlet, tbls, mappings, opts, 'useVirtual', useVirtual);

    matrices = output.matrices;
    bcvals   = output.bcvals;
    extra    = output.extra;
    
    % We enforce the Dirichlet boundary conditions as Lagrange multipliers
    if ~isempty(bcstruct.bcneumann)
        error('not yet implemented');
    else
        nf_num = nodefacetbl.num;
        extflux = zeros(nf_num, 1);
    end

    if isempty(src)
        src = zeros(celltbl.num, 1);
    end
    
    
    fullrhs{1} = extflux;
    fullrhs{2} = src;
    fullrhs{3} = bcvals;
    
    if (opt.onlyAssemblyMatrices | opt.addAssemblyMatrices | opt.addAdOperators)
        
        matrices.fullrhs = fullrhs;
        assembly.matrices = matrices;
        assembly.nKg = extra.nKg;
        
        if opt.onlyAssemblyMatrices
            return
        end
        
    end
    
    % We reduced the system (shur complement) using invA11
    % We obtain system of the form
    %
    % B*u = rhs
    %
    % where
    %
    % B = [[B11, B12];
    %      [B21, B22]];
    %
    % u = [p (pressure at celltbl);
    %      lagmult];
    %
    % rhs = [-A21*invA11*extflux;  +  [src;
    %        -D'*invA11*extflux  ]     bcvals]
    
    invA11 = matrices.invA11; 
    A12 = matrices.A12;
    A21 = matrices.A21;
    A22 = matrices.A22;
    D = matrices.D;
    
    B11 = A22 - A21*invA11*A12;
    B12 = A21*invA11*D;
    B21 = -D'*invA11*A12;
    B22 = D'*invA11*D;


    B = [[B11, B12]; ...
         [B21, B22]];
    
    adrhs{1} = -A21*invA11*extflux + src; 
    adrhs{2} = -D'*invA11*extflux + bcvals;
    
    rhs = vertcat(adrhs{:});
    
    assembly.B = B;
    assembly.rhs = rhs;
    
    if opt.addAdOperators
        
        adB = cell(2, 2);
        adB{1, 1} = B11;
        adB{2, 1} = B21;
        adB{1, 2} = B12;
        adB{2, 2} = B22;
        
        adoperators.B     = adB;
        adoperators.rhs   = adrhs;        
        
        % Setup fluid flux operator
        mpfaKgrad = setupMpfaFlux(G, assembly, tbls);
        fluxop = @(p) fluxFunc(p, mpfaKgrad);
        
        adoperators.fluxop = fluxop;
        
        assembly.adoperators = adoperators;
        
    end
    
end

function flux = fluxFunc(p, mpfaKgrad)
   flux = mpfaKgrad*p;
end
