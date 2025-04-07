function C = setupStiffnessTensor(prop, tbls, mappings, varargin)
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


    opt = struct('useVirtual', false);
    opt = merge_options(opt, varargin{:});

    useVirtual = opt.useVirtual;
    
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
    
    % constructiontypes = {'using_change_of_basis', ...
    %                      'general_voigt_construction'};
    constructiontype = 'using_change_of_basis';
    
    switch constructiontype

      case 'using_change_of_basis'

        mu_     = unique(mu);
        lambda_ = unique(lambda);

        if (numel(mu_) == 1 && numel(lambda_) == 1)
            isscalar = true;
            mu     = mu_;
            lambda = lambda_;
        else
            isscalar = false;
        end

        switch dim
            
          case 2
            
            M = [1 0 0 0;
                 0 0 1 1;
                 0 0 1 -1;
                 0 1 0 0];
            
          case 3

            M = [[1, 0, 0, 0, 0, 0, 0 , 0 , 0];
                 [0, 0, 0, 1, 0, 0, 1 , 0 , 0];
                 [0, 0, 0, 0, 1, 0, 0 , 1 , 0];
                 [0, 0, 0, 1, 0, 0, -1, 0 , 0];
                 [0, 1, 0, 0, 0, 0, 0 , 0 , 0];
                 [0, 0, 0, 0, 0, 1, 0 , 0 , 1];
                 [0, 0, 0, 0, 1, 0, 0 , -1, 0];
                 [0, 0, 0, 0, 0, 1, 0 , 0 , -1];
                 [0, 0, 1, 0, 0, 0, 0 , 0 , 0]];


        end

        invM = inv(M);

        if isscalar
            
            switch dim
                
              case 2
                
                C = [(lambda + 2*mu), lambda         , 0   , 0;
                     lambda         , (lambda + 2*mu), 0   , 0;
                     0              , 0              , 2*mu, 0;
                     0              , 0              , 0   , 2*mu];
                
              case 3

                C = [[(lambda + 2*mu), lambda         , lambda         , 0   , 0   , 0   , 0   , 0   , 0];
                     [lambda         , (lambda + 2*mu), lambda         , 0   , 0   , 0   , 0   , 0   , 0];
                     [lambda         , lambda         , (lambda + 2*mu), 0   , 0   , 0   , 0   , 0   , 0];
                     [0              , 0              , 0              , 2*mu, 0   , 0   , 0   , 0   , 0];
                     [0              , 0              , 0              , 0   , 2*mu, 0   , 0   , 0   , 0];
                     [0              , 0              , 0              , 0   , 0   , 2*mu, 0   , 0   , 0];
                     [0              , 0              , 0              , 0   , 0   , 0   , 2*mu, 0   , 0];
                     [0              , 0              , 0              , 0   , 0   , 0   , 0   , 2*mu, 0];
                     [0              , 0              , 0              , 0   , 0   , 0   , 0   , 0   , 2*mu]];
            end

            C = M*C*invM;
            C = reshape(C', [], 1);
            
            map = TensorMap();
            map.fromTbl  = vec1212tbl;
            map.toTbl    = cellvec1212tbl;
            map.mergefds = {'vec11', 'vec12', 'vec21', 'vec22'};

            if useVirtual

                map.pivottbl = cellvec1212tbl;

                [vec, i] = ind2sub([vec1212tbl.num, celltbl.num], (1 : cellvec1212tbl.num)');
                map.dispind1 = vec;
                map.dispind2 = (1 : cellvec1212tbl.num)';
                map.issetup = true;
                
            else
                
                map = map.setup();
                
            end
            
            C = map.eval(C);
            
        else
            
            prod = TensorProd();
            prod.tbl1 = vec1212tbl;
            prod.tbl2 = celltbl;
            prod.tbl3 = cellvec1212tbl;

            if useVirtual

                prod.pivottbl = cellvec1212tbl;
                [vec, i] = ind2sub([vec1212tbl.num, celltbl.num], (1 : cellvec1212tbl.num)');                
                prod.dispind1 = vec;
                prod.dispind2 = i;
                prod.dispind3 = (1 : cellvec1212tbl.num)';
                prod.issetup = true;
                
            else
                
                prod = prod.setup();
                
            end

            u = zeros(vec1212tbl.num, 1);
            u(vec1212tbl.get('vec11') == vec1212tbl.get('vec21') & vec1212tbl.get('vec12') == vec1212tbl.get('vec22') )= 1;

            muC = prod.eval(u, 2*mu);

            u = zeros(vec1212tbl.num, 1);
            u(vec1212tbl.get('vec11') == 1 & vec1212tbl.get('vec21') == 1) = 1;

            lambdaC = prod.eval(u, lambda);

            tC = lambdaC + muC;
            
            tconv = TensorConvert;
            tconv.fromTbl           = vec12tbl;
            tconv.toTbl             = vec12tbl;
            tconv.pivottbl          = vec1212tbl;
            tconv.replacefdsFromTbl = {{'vec1', 'vec21'}, {'vec2', 'vec22'}};
            tconv.replacefdsToTbl   = {{'vec1', 'vec11'}, {'vec2', 'vec12'}};
            tconv = tconv.setup();

            Mv    = tconv.convert(M);
            invMv = tconv.convert(invM);

            prod = TensorProd();
            prod.tbl1        = vec1212tbl;
            prod.tbl2        = cellvec1212tbl;
            prod.tbl3        = cellvec1212tbl;
            prod.replacefds1 = {{'vec11', 'redvec1'}, {'vec12', 'redvec2'}};
            prod.replacefds2 = {{'vec21', 'redvec1'}, {'vec22', 'redvec2'}};
            prod.reducefds   = {'redvec1', 'redvec2'};

            if useVirtual

                gen = CrossIndexArrayGenerator();
                gen.tbl1        = vec1212tbl;
                gen.tbl2        = vec12tbl;
                gen.replacefds2 = {{'vec1', 'vec31'}, {'vec2', 'vec32'}};
                gen.mergefds    = {};
                gen.opts        = {'optpureproduct', true, 'virtual', useVirtual};
                
                vec121212tbl = gen.eval();

                cellvec121212tbl = crossIndexArray(celltbl, vec121212tbl, {}, 'optpureproduct', true, 'virtual', useVirtual);
                
                prod.pivottbl = cellvec121212tbl;

                [vec3, vec2, vec1, i] = ind2sub([vec12tbl.num, vec12tbl.num, vec12tbl.num, celltbl.num], (1 : cellvec121212tbl.num)');
                prod.dispind1 = sub2ind([vec12tbl.num, vec12tbl.num], vec1, vec2);
                prod.dispind2 = sub2ind([vec12tbl.num, vec12tbl.num, celltbl.num], vec2, vec3, i);;
                prod.dispind3 = sub2ind([vec12tbl.num, vec12tbl.num, celltbl.num], vec1, vec3, i);;
                prod.issetup = true;
                
            else
                
                prod = prod.setup();
                
            end
            
            C = prod.eval(invMv, tC);

            prod = TensorProd();
            prod.tbl1        = vec1212tbl;
            prod.tbl2        = cellvec1212tbl;
            prod.tbl3        = cellvec1212tbl;
            prod.replacefds1 = {{'vec21', 'redvec1'}, {'vec22', 'redvec2'}};
            prod.replacefds2 = {{'vec11', 'redvec1'}, {'vec12', 'redvec2'}};
            prod.reducefds   = {'redvec1', 'redvec2'};

            if useVirtual
                
                prod.pivottbl = cellvec121212tbl;

                [vec3, vec2, vec1 i] = ind2sub([vec12tbl.num, vec12tbl.num, vec12tbl.num, celltbl.num], (1 : cellvec121212tbl.num)');
                prod.dispind1 = sub2ind([vec12tbl.num, vec12tbl.num], vec2, vec3);
                prod.dispind2 = sub2ind([vec12tbl.num, vec12tbl.num, celltbl.num], vec1, vec2, i);;
                prod.dispind3 = sub2ind([vec12tbl.num, vec12tbl.num, celltbl.num], vec1, vec3, i);;
                prod.issetup = true;
                
            else
                
                prod = prod.setup();
                
            end

            C = prod.eval(Mv, C);
            
        end
        
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
