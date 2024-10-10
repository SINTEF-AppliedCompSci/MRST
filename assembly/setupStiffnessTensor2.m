function C = setupStiffnessTensor2(prop, tbls)
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


    celltbl        =  tbls.celltbl;
    vectbl         =  tbls.vectbl;
    vec12tbl       =  tbls.vec12tbl;
    cellvec1212tbl =  tbls.cellvec1212tbl;
    vec1212tbl     =  tbls.vec1212tbl;
    
    dim = vectbl.num;
    
    lambda = prop.lambda;
    mu = prop.mu;
    
    vdim = dim*(dim + 1)/2;
    avdim = dim*dim - vdim;
    
    constructiontypes = {'direct_lambda_mu_construction', ...
                        'general_voigt_construction'};
    % constructiontype = 'general_voigt_construction';
    % constructiontype = 'direct_lambda_mu_construction';
    constructiontype = 'using_change_of_basis';
    
    switch constructiontype

      case 'using_change_of_basis'

        assert(dim == 2, 'only for 2d at the moment!');
        
        mu     = unique(mu);
        lambda = unique(lambda);

        M = [1 0 0 0;
             0 0 1 1;
             0 0 1 -1;
             0 1 0 0];

        C = [(lambda + 2*mu), lambda         , 0   , 0;
             lambda         , (lambda + 2*mu), 0   , 0;
             0              , 0              , 2*mu, 0;
             0              , 0              , 0   , 2*mu];

        C = M*C*inv(M);

        C = reshape(C', [], 1);
        
        map = TensorMap();
        map.fromTbl = vec1212tbl;
        map.toTbl = cellvec1212tbl;
        map.mergefds = {'vec11', 'vec12', 'vec21', 'vec22'};
        map = map.setup();

        C = map.eval(C);
        
      case 'direct_lambda_mu_construction'
        
        % see formula https://en.wikipedia.org/wiki/Hooke%27s_law
        clear mutbl
        mutbl.vec11 = vec12tbl.get('vec1');
        mutbl.vec21 = vec12tbl.get('vec1');
        mutbl.vec12 = vec12tbl.get('vec2');
        mutbl.vec22 = vec12tbl.get('vec2');
        mutbl1 = IndexArray(mutbl);

        clear mutbl
        mutbl.vec11 = vec12tbl.get('vec1');
        mutbl.vec22 = vec12tbl.get('vec1');
        mutbl.vec12 = vec12tbl.get('vec2');
        mutbl.vec21 = vec12tbl.get('vec2');
        mutbl2 = IndexArray(mutbl);

        map = TensorMap();
        map.fromTbl  = mutbl1;
        map.toTbl    = vec1212tbl;
        map.mergefds = {'vec11', 'vec12', 'vec21', 'vec22'};
        map = map.setup();

        mucoef1 = map.eval(ones(mutbl1.num, 1));

        map = TensorMap();
        map.fromTbl  = mutbl2;
        map.toTbl    = vec1212tbl;
        map.mergefds = {'vec11', 'vec12', 'vec21', 'vec22'};
        map = map.setup();

        mucoef2 = map.eval(ones(mutbl2.num, 1));

        mucoef = (mucoef1 + mucoef2);

        prod = TensorProd();
        prod.tbl1 = celltbl;
        prod.tbl2 = vec1212tbl;
        prod.tbl3 = cellvec1212tbl;
        prod = prod.setup();

        Cmu = prod.eval(mu, mucoef);

        diagtbl.vec1 = (1 : dim)';
        diagtbl.vec2 = (1 : dim)';
        diagtbl = IndexArray(diagtbl);

        gen = CrossIndexArrayGenerator;
        gen.tbl1 = diagtbl;
        gen.tbl2 = diagtbl;
        gen.replacefds1 = {{'vec1', 'vec11'}, {'vec2', 'vec12'}};
        gen.replacefds2 = {{'vec1', 'vec21'}, {'vec2', 'vec22'}};

        lambdatbl = gen.eval();

        celllambdatbl = crossIndexArray(celltbl, lambdatbl, {});

        map = TensorMap();
        map.fromTbl = celltbl;
        map.toTbl = celllambdatbl;
        map.mergefds = {'cells'};
        map = map.setup();

        lambda = map.eval(lambda);

        map = TensorMap();
        map.fromTbl = celllambdatbl;
        map.toTbl = cellvec1212tbl;
        fds = {'vec11', 'vec12', 'vec21', 'vec22', 'cells'};
        map.mergefds = fds;
        map = map.setup();

        Clambda = map.eval(lambda);

        C = Cmu + Clambda;
        
      case 'general_voigt_construction'

        error('not fully functional yet');
        
        Cvoigt = mu*eye(vdim);
        Z1 = zeros(dim, vdim - dim);
        Z2 = zeros(vdim - dim);
        Cvoigt = [[lambda*ones(dim), Z1]; ...
                  [Z1'             , Z2]] ...
                 + Cvoigt;

        % Asymmetric part
        Casym = mu*eye(avdim);

        n1 = size(Cvoigt, 1);
        n2 = size(Casym, 1);

        Z = zeros(n1, n2);

        C = [[Cvoigt, Z];...
             [Z', Casym];
            ];

        C = reshape(C', [], 1);

        % We dispatch C
        map = TensorMap();
        map.fromTbl = vec1212tbl;
        map.toTbl = cellvec1212tbl;
        fds = {'vec11', 'vec12', 'vec21', 'vec22'};
        map.mergefds = fds;
        map = map.setup();
        
        C = map.eval(C);
        
        % Change of basis :  
        % mapping A -> A + A' and A -> A - A'
        % follows indexing of vec22tbl

        voigttbl.voigt = (1 : vdim)';
        voigttbl.num = numel(voigttbl.voigt);

        avoigttbl.avoigt = (1 : avdim)';
        avoigttbl.num = numel(avoigttbl.avoigt);

        vec12voigttbl.coldim = vec12tbl.coldim;
        vec12voigttbl.rowdim = vec12tbl.rowdim;
        switch dim
          case 2
            vec12voigttbl.voigt  = [1; 3; 3; 2];
          case 3
            vec12voigttbl.voigt  = [1; 6; 5; 6; 2; 4; 5; 4; 3];
        end

        vec12voigttbl.num = numel(vec12voigttbl.coldim);

        % to index the avoigt we the same ordering as voigt, just skipping the diagonal
        switch dim
          case 2
            vec12avoigttbl.avoigt = [1; 1];
            vec12avoigttbl.coldim = [1; 2];
            vec12avoigttbl.rowdim = [2; 1];
          case 3
            vec12avoigttbl.avoigt = [1; 2; 3; 1; 2; 3];
            vec12avoigttbl.coldim = [2; 1; 1; 3; 3; 2];
            vec12avoigttbl.rowdim = [3; 3; 2; 2; 1; 1];
        end

        vec12avoigttbl.num = numel(vec12avoigttbl.avoigt);

        prod = TensorProd();
        prod.tbl1 = vec12voigttbl;
        prod.tbl2 = voigttbl;
        prod.tbl3 = vec12tbl;
        prod.reducefds = {'voigt'};
        prod = prod.setup();

        V_T = SparseTensor();
        V_T = V_T.setFromTensorProd(ones(vec12voigttbl.num, 1), prod);
        V = V_T.getMatrix();

        prod = TensorProd();
        prod.tbl1 = vec12avoigttbl;
        prod.tbl2 = avoigttbl;
        prod.tbl3 = vec12tbl;
        prod.reducefds = {'avoigt'};
        prod = prod.setup();

        AV_T = SparseTensor();
        coef = [ones(avdim, 1); -ones(avdim, 1)]; 
        AV_T = AV_T.setFromTensorProd(coef, prod);
        AV = AV_T.getMatrix();

        M1 = [V, AV]';
        M1 = full(M1);

        % We add an extra multiplication for the coef 2 on diagonal for Voigt.
        d = [ones(dim, 1); 0.5*ones(dim*dim - dim, 1)];
        Vc = diag(d);

        M2 = M1'*Vc;

        % we write M1 and M2 as elements of col2row2tbl
        M1 = reshape(M1', [], 1);
        M2 = reshape(M2', [], 1);

        prod = TensorProd();
        prod.tbl1 = cellcol2row2tbl;
        prod.tbl2 = col2row2tbl;
        prod.tbl3 = cellcol2row2tbl;
        prod.replacefds1 = {{'coldim2', 'coldim'}, {'rowdim2', 'rowdim'}};
        prod.replacefds2 = {{'coldim1', 'coldim'}, {'rowdim1', 'rowdim'}};
        prod.mergefds = {'coldim', 'rowdim'};

        prod = prod.setup();

        CM1 = prod.eval(C, M1);
        
        prod = TensorProd();
        prod.tbl1 = col2row2tbl;
        prod.tbl2 = cellcol2row2tbl;
        prod.tbl3 = cellcol2row2tbl;
        prod.replacefds1 = {{'coldim2', 'coldim'}, {'rowdim2', 'rowdim'}};
        prod.replacefds2 = {{'coldim1', 'coldim'}, {'rowdim1', 'rowdim'}};
        prod.mergefds = {'coldim', 'rowdim'};

        prod = prod.setup();

        M2CM1 = prod.eval(M2, CM1);
        C = M2CM1;
        
        dotest = false;
        if dotest
            prod = TensorProd();
            prod.tbl1 = col2row2tbl;
            prod.tbl2 = vec12tbl;
            prod.replacefds1 = {{'coldim1', 'coldim'}, {'rowdim1', 'rowdim'}};
            prod.replacefds2 = {{'coldim', 'coldim2'}, {'rowdim', 'rowdim2'}};
            prod.reducefds = {'coldim2', 'rowdim2'};
            prod.tbl3 = vec12tbl;
            prod = prod.setup();

            C_T = SparseTensor('matlabsparse', true);
            C_T = C_T.setFromTensorProd(C, prod);
            
            printTensor(C_T); 
        end 
        
      otherwise
        error('constructiontype not recognized');
    end

end
