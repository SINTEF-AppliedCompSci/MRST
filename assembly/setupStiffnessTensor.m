function C = setupStiffnessTensor(prop, tbls)
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


    celltbl         =  tbls.celltbl;
    coltbl          =  tbls.coltbl;
    colrowtbl       =  tbls.colrowtbl;
    cellcol2row2tbl =  tbls.cellcol2row2tbl;
    col2row2tbl     =  tbls.col2row2tbl;
    
    dim = coltbl.num;
    
    lambda = prop.lambda;
    mu = prop.mu;
    
    vdim = dim*(dim + 1)/2;
    avdim = dim*dim - vdim;

    
    constructiontypes = {'direct_lambda_mu_construction', ...
                        'general_voigt_construction'};
    % constructiontype = 'general_voigt_construction';
    constructiontype = 'direct_lambda_mu_construction';

    switch constructiontype
      case 'direct_lambda_mu_construction'

        mutbl.coldim1 = colrowtbl.get('coldim');
        mutbl.coldim2 = colrowtbl.get('coldim');
        mutbl.rowdim1 = colrowtbl.get('rowdim');
        mutbl.rowdim2 = colrowtbl.get('rowdim');
        mutbl = IndexArray(mutbl);

        cellmutbl = crossIndexArray(celltbl, mutbl, {});

        map = TensorMap();
        map.fromTbl = celltbl;
        map.toTbl = cellmutbl;
        map.mergefds = {'cells'};
        map = map.setup();

        Cmu = map.eval(2*mu);

        map = TensorMap();
        map.fromTbl = cellmutbl;
        map.toTbl = cellcol2row2tbl;
        fds = {'coldim1', 'coldim2', 'rowdim1', 'rowdim2', 'cells'};
        map.mergefds = fds;
        map = map.setup();

        Cmu = map.eval(Cmu);

        diagtbl.coldim = (1 : dim)';
        diagtbl.rowdim = (1 : dim)';
        diagtbl = IndexArray(diagtbl);

        fds = {{'rowdim', {'rowdim1', 'rowdim2'}}, ...
               {'coldim', {'coldim1', 'coldim2'}}};
        lambdatbl = crossIndexArray(diagtbl, diagtbl, {}, 'crossextend', fds);

        celllambdatbl = crossIndexArray(celltbl, lambdatbl, {});

        map = TensorMap();
        map.fromTbl = celltbl;
        map.toTbl = celllambdatbl;
        map.mergefds = {'cells'};
        map = map.setup();

        lambda = map.eval(lambda);

        map = TensorMap();
        map.fromTbl = celllambdatbl;
        map.toTbl = cellcol2row2tbl;
        fds = {'coldim1', 'coldim2', 'rowdim1', 'rowdim2', 'cells'};
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
        map.fromTbl = col2row2tbl;
        map.toTbl = cellcol2row2tbl;
        fds = {'coldim1', 'coldim2', 'rowdim1', 'rowdim2'};
        map.mergefds = fds;
        map = map.setup();
        
        C = map.eval(C);
        
        % Change of basis :  
        % mapping A -> A + A' and A -> A - A'
        % follows indexing of colrowtbl

        voigttbl.voigt = (1 : vdim)';
        voigttbl.num = numel(voigttbl.voigt);

        avoigttbl.avoigt = (1 : avdim)';
        avoigttbl.num = numel(avoigttbl.avoigt);

        colrowvoigttbl.coldim = colrowtbl.coldim;
        colrowvoigttbl.rowdim = colrowtbl.rowdim;
        switch dim
          case 2
            colrowvoigttbl.voigt  = [1; 3; 3; 2];
          case 3
            colrowvoigttbl.voigt  = [1; 6; 5; 6; 2; 4; 5; 4; 3];
        end

        colrowvoigttbl.num = numel(colrowvoigttbl.coldim);

        % to index the avoigt we the same ordering as voigt, just skipping the diagonal
        switch dim
          case 2
            colrowavoigttbl.avoigt = [1; 1];
            colrowavoigttbl.coldim = [1; 2];
            colrowavoigttbl.rowdim = [2; 1];
          case 3
            colrowavoigttbl.avoigt = [1; 2; 3; 1; 2; 3];
            colrowavoigttbl.coldim = [2; 1; 1; 3; 3; 2];
            colrowavoigttbl.rowdim = [3; 3; 2; 2; 1; 1];
        end

        colrowavoigttbl.num = numel(colrowavoigttbl.avoigt);

        prod = TensorProd();
        prod.tbl1 = colrowvoigttbl;
        prod.tbl2 = voigttbl;
        prod.tbl3 = colrowtbl;
        prod.reducefds = {'voigt'};
        prod = prod.setup();

        V_T = SparseTensor();
        V_T = V_T.setFromTensorProd(ones(colrowvoigttbl.num, 1), prod);
        V = V_T.getMatrix();

        prod = TensorProd();
        prod.tbl1 = colrowavoigttbl;
        prod.tbl2 = avoigttbl;
        prod.tbl3 = colrowtbl;
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
            prod.tbl2 = colrowtbl;
            prod.replacefds1 = {{'coldim1', 'coldim'}, {'rowdim1', 'rowdim'}};
            prod.replacefds2 = {{'coldim', 'coldim2'}, {'rowdim', 'rowdim2'}};
            prod.reducefds = {'coldim2', 'rowdim2'};
            prod.tbl3 = colrowtbl;
            prod = prod.setup();

            C_T = SparseTensor('matlabsparse', true);
            C_T = C_T.setFromTensorProd(C, prod);
            
            printTensor(C_T); 
        end 
        
      otherwise
        error('constructiontype not recognized');
    end

end
