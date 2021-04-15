classdef SparseTensor
%
% SYNOPSIS:
%   SparseTensor(varargin)
%
% DESCRIPTION:
%
% Class for Sparse Tensor, which can be seen as a sparse matrix. However, using
% our own implementation enables us to use AD variables for the coefficients
% (stored in tensval) while the matlab format only accepts scalar values.
%
% We also keep track of the spaces the tensor maps from and to (they are given
% as IndexArray). The tensor maps from fromTbl to toTbl
%
% We overload the multipication and addition operations.
%
% An essential method is setFromTensorProd which sets a linear mapping from a
% TensorProd (denoted prod) and a vector (denoted vec). Indeed, the mapping u ->
% v given by 
%
%    v = prod.eval(vec, u) 
%    
% is linear. The values of u follow the indexing given by the IndexArray fromTbl
% and the values of v follow the indexing given by the IndexArray toTbl.
%
% Any linear mapping can be written in this form, for some (non-unique) prod and
% vec. In the assembly of discretization methods, the linear mappings are
% typically presented in this form.
%
% PARAMETERS:
%   varargin - 
%
% RETURNS:
%   class instance
%
% EXAMPLE:
%   `computeVagTrans`
%
% SEE ALSO:
%   `IndexArray`, `TensorProd`.
    
    properties
    
        tensval      % Value of the coefficients of the tensor
        col          % Column index for tensval (indexing corresponds to local indexing of fromTbl)
        row          % Row index for tensval (indexing corresponds to local indexing of toTbl)
        fromTbl      % Argument space
        toTbl        % Result space
        matrix       % Sparse matrix in Matlab format
        matlabsparse % If set to true the matlab sparse matrix format will be set and
                     % used by the multiplication and addition
    end
   
    methods
        
        function tensor = SparseTensor(varargin)
            opt = struct('prod'    , [], ...
                         'value'   , [], ...
                         'matlabsparse', false, ...
                         'argindex', 1);
            opt = merge_options(opt, varargin{:}); 
            
            if ~isempty(opt.prod)
                tensor = tensor.setFromTensorProd(opt.value, opt.prod, ...
                                                             'argindex', ...
                                                             opt.argindex);
            end
            
            tensor.matrix = [];
            tensor.matlabsparse = opt.matlabsparse;
            
        end

        function tensor = setMatrix(tensor)
        % set the sparse matrix in Matlab format
            val     = tensor.tensval;  
            col     = tensor.col;
            row     = tensor.row;
            fromTbl = tensor.fromTbl; 
            toTbl   = tensor.toTbl;   
            
            if isa(val, 'ADI')
                val = val.value;
            end
            
            mat = sparse(row, col, val, toTbl.num, fromTbl.num);
            
            tensor.matrix = mat;
            
        end
        
        function mat = getMatrix(tensor)
        % get the sparse matrix in Matlab format
            if  isempty(tensor.matrix) 
                tensor = tensor.setMatrix();
            end
            mat = tensor.matrix;

        end
        
        function tensor = setFromTensorProd(tensor, vec, prod, varargin)
            opts = struct('argindex', 1);            
            opts = merge_options(opts, varargin{:}); 
            argindex = opts.argindex;
            
            if ~(prod.issetup)
                prod = prod.setup();
            end
            
            dispind1 = prod.dispind1;
            dispind2 = prod.dispind2;
            dispind3 = prod.dispind3;
            
            switch argindex
              case 1
                fromtbl = prod.tbl2;
                dispind = dispind2;
                vec   = vec(dispind1);
              case 2
                fromtbl = prod.tbl1;
                dispind = dispind1;
                vec   = vec(dispind2);
            end
            
            redind  = dispind3;
            fromTbl = fromtbl;
            toTbl   = prod.tbl3;

            tensor.tensval = vec;
            tensor.col     = dispind;
            tensor.row     = redind;
            tensor.fromTbl = fromTbl;
            tensor.toTbl   = toTbl;
            
            if tensor.matlabsparse
                tensor = tensor.setMatrix();
            end
            
        end

        function tensor = setFromTensorMap(tensor, map)

            if ~(map.issetup)
                map = map.setup();
            end
            
            tensor.tensval = ones(size(map.dispind1, 1), 1);
            tensor.col     = map.dispind1;
            tensor.row     = map.dispind2;
            tensor.fromTbl = map.fromTbl;
            tensor.toTbl   = map.toTbl;
            
            if tensor.matlabsparse
                tensor = tensor.setMatrix();
            end
            
        end

        function ttens = transpose(tens)
                
            ttens = SparseTensor();
            ttens.tensval     = tens.tensval;
            ttens.col     = tens.row;
            ttens.row     = tens.col;
            ttens.fromTbl = tens.toTbl;
            ttens.toTbl   = tens.fromTbl;
            
        end
        
        function tens = mtimes(tens1, tens2)
            
            if isa(tens1, 'double')
                c    = tens1;
                tens = tens2;
                
                if  ~isempty(tens.matrix)
                    tens.matrix = c*tens.matrix;
                    if ~tens.matlabsparse
                        tens.tensval = c*tens.tensval;
                    end
                    return
                end
                tens.tensval = c*tens.tensval;
                
            else
                
                tens = SparseTensor();
                tens.fromTbl = tens2.fromTbl;
                tens.toTbl   = tens1.toTbl;

                if ~isempty(tens1.matrix) & ~isempty(tens2.matrix)
                    mat1 = tens1.matrix;
                    mat2 = tens2.matrix;
                    mat = mat1*mat2;
                    tens.matrix = mat;
                    if ~(tens1.matlabsparse) | ~(tens2.matlabsparse)
                        [row, col, tensval] = find(mat);
                        tens.col = col;                 
                        tens.row = row;
                        tens.tensval = tensval;
                    end
                    return
                end
                
                row1 = tens1.row;
                col1 = tens1.col;
                tensval1 = tens1.tensval;
                
                row2 = tens2.row;
                col2 = tens2.col;
                tensval2 = tens2.tensval;
                
                tbl1.A = row1;
                tbl1.B = col1;
                tbl1 = IndexArray(tbl1);
                
                tbl2.B = row2;
                tbl2.C = col2;
                tbl2 = IndexArray(tbl2);
                
                [tbl3, indstruct] = crossIndexArray(tbl1, tbl2, {'B'});
                
                tensval1 = tbldispatch1(tensval1, indstruct);
                tensval2 = tbldispatch2(tensval2, indstruct);
                
                tensval = tensval1.*tensval2;
                
                fds = {'A', 'C'};
                tbl = projIndexArray(tbl3, fds);
                
                map = TensorMap();
                map.fromTbl = tbl3;
                map.toTbl = tbl;
                map.mergefds = fds;
                map = map.setup();
                
                tensval = map.eval(tensval);
                
                
                tens.tensval = tensval;
                tens.col = tbl.get('C');
                tens.row = tbl.get('A');
            end
            
        end
        
        function tens = plus(tens1, tens2)
            
            tens = SparseTensor();
            tens.fromTbl = tens1.fromTbl; % = tens2.fromTbl
            tens.toTbl   = tens1.toTbl; % = tens2.toTbl

            if ~isempty(tens1.matrix) & ~isempty(tens2.matrix)
                mat1 = tens1.matrix;
                mat2 = tens2.matrix;
                mat = mat1 + mat2;
                tens.matrix = mat;
                if ~(tens1.matlabsparse) | ~(tens2.matlabsparse)
                    [row, col, tensval] = find(mat);
                    tens.col = col;                 
                    tens.row = row;
                    tens.tensval = tensval;
                end
                return
            end
            
            row1 = tens1.row;
            col1 = tens1.col;
            tensval1 = tens1.tensval;
            
            row2 = tens2.row;
            col2 = tens2.col;
            tensval2 = tens2.tensval;
            
                        
            tbl1.A = row1;
            tbl1.B = col1;
            tbl1.num = numel(tbl1.A);
            
            tbl2.A = row2;
            tbl2.B = col2;
            tbl2.num = numel(tbl2.B);
            
            rowcol = [[row1; row2], [col1; col2]];
            rowcol = unique(rowcol, 'rows');
            
            tbl.A = rowcol(:, 1);
            tbl.B = rowcol(:, 2);
            tbl.num = size(rowcol, 1);

            map = TensorMap();
            map.fromTbl = tbl1;
            map.toTbl = tbl;
            map.mergefds = {'A', 'B'};
            map = map.setup();
            
            tensval1 = map.eval(tensval1);
            
            map = TensorMap();
            map.fromTbl = tbl2;
            map.toTbl = tbl;
            map.mergefds = {'A', 'B'};
            map = map.setup();
            
            tensval2 = map.eval(tensval2);
            
            tensval = tensval1 + tensval2;
            
            tens.tensval = tensval;
            tens.row = tbl.A;
            tens.col = tbl.B;
            
        end
        

    end
    
        
   
end

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
