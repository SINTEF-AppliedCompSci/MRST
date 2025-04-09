classdef TensorProd
%
%
% SYNOPSIS:
%   TensorProd(varargin)
%
% DESCRIPTION:
%
% Class used to define a tensor product between two spaces. The two spaces are
% described using IndexArray. The product is defined by the fields that are
% going to be merged and reduced. We use this class as follows: First, we
% instantiate a product (prod = TensorProd), then we set the mergefds and
% reducefds index names (for example, prod.mergefds = {'A'}, prod.reducefds =
% {'B', 'C'} for example). Finally, we setup the product by running prod =
% prod.setup(). In this last stage, an internal representation of the product
% (using local linear indices) is created and the product of two vectors (w =
% prod.eval(u, v)) will be computed from that efficiently. This representation
% is given by the index vectors dispind1 dispind2, dispind3 and the IndexArray
% pivottbl. In particular, the sparsity for an efficient computation of the
% product is computed (see paper on Tensor Product assembly).
%
% We explain the meaning of merging and reducing along the indices given by
% mergefds and reducefds. Let us write the multi-indexed vectors:
%
%  u[i, j, k] in IndexArray tbl1,
%  v[i, j, l] in IndexArray tbl2 
%   
% where :
%  * the first index i corresponds to the mergefds indices (i is index in 'A' for example),
%  * the second index j corresponds to the reducefds indices (j is a (multi)-index in {'A', 'B'} for example).
%  
% Then the result w of the product of u and v, that is setup here, is given by
%   
%    w[i, k, l] = sum_{i, j} (u[i, j, k] v[i, j, l])
%
% The index i is merged (it appears in u, v and also in w). The index j is
% reduced (it appears in u, v but is summed up and does not appear in w).
%
% In the setup of the product (before calling prod = prod.setup()), the
% fieldnames of tbl1 and tbl2 can be changed using prod.replacefds1 and
% prod.replacefds2. This option is used very often as different names for the
% same indexing appear naturally.
%
% If the pivot space (pivottbl property) is given before the setup, then it is possible to dispatch in a index array
% tbl3 which does not corresponds to the otherwise detected sparsity of the product. It is implemented for convenience,
% to avoid to have to create an extra TensorMap to do this this operation.
%
% RETURNS:
%   class instance
%
% EXAMPLE:
%   `computeVagTrans`
%
% SEE ALSO:
%   `IndexArray`, `SparseTensor`.
    
    properties
        tbl1         % Table for first argument
        tbl2         % Table for second argument
        tbl3         % Table for result 
        reducefds    % Field names for reduction
        mergefds     % Field names for merging
        
        replacefds1  % Possibility to change field names of first table
                     % (before setting up product)
        replacefds2  % Possibility to change field names of second table
                     % (before setting up product)
        replacefds3  % Possibility to change field names of product table (after setting
                     % up product) for the default tbl3 (not applied if tbl3 is
                     % given as input). This option is little used, as it can be
                     % avoid by using only replacefds1 and replacefds2

        reducefds1   % Field names that belongs to tbl1 and not tbl2 that
                     % will be reduced
        reducefds2   % Field names that belongs to tbl2 and not tbl1 that
                     % will be reduced
        
        % Structures used for the computation. Initialized only once.
        pivottbl   % Structure for the expanded product
        
        % Dispatch indices 
        dispind1 % from tbl1 to pivottbl 
        dispind2 % from tbl2 to pivottbl
        dispind3 % from tbl3 to pivottbl
        
        issetup      % Flag is set to true is product has been set up.
        settbl3      % Flag is set to true if the resulting product table is created
        setpivottbl3 % Flag is set to true if the pivottbl table is created

        chunksize = 500000; % Chunk size for the computation of the product

        needIntersectTbl3 = false; % For the automatically generated pivot space, the index array tbl3 must be a projection of the pivot
                                   % space. For some cases, one may want to use a tbl3 which is not a projection
                                   % of the pivot space, and will require an extra mapping. Then, this option should be
                                   % set to true, so that the intersection is computed.

        useMex = true;
        
    end
   
    methods
        
        function prod = TensorProd(varargin)
            
            opts = struct('tbl1'       , [], ...
                          'tbl2'       , [], ...
                          'tbl3'       , [], ...
                          'replacefds1', [], ...
                          'replacefds2', [], ...
                          'replacefds3', [], ...
                          'reducefds'  , [], ...
                          'reducefds1' , [], ...
                          'reducefds2' , [], ...
                          'mergefds'   , []);
            
            prod = merge_options(prod, varargin{:}); 
            
            if isempty(prod.reducefds)
                prod.reducefds = {};
            end
            
            if isempty(prod.reducefds1)
                prod.reducefds1 = {};
            end
            
            if isempty(prod.reducefds2)
                prod.reducefds2 = {};
            end

            if isempty(prod.mergefds)
                prod.mergefds = {};
            end            
            
            if isempty(prod.tbl3)
                prod.settbl3 = true;
            else
                prod.settbl3 = false;
            end
            
            prod.pivottbl = [];
            prod.issetup = false;
            
        end
        
        function prod = setup(prod)
            
            reducefds  = prod.reducefds;
            reducefds1 = prod.reducefds1;
            reducefds2 = prod.reducefds2;
            mergefds   = prod.mergefds;

            crossfds  = {reducefds{:}, mergefds{:}};
            
            tbl1 = prod.tbl1;
            tbl2 = prod.tbl2;
            
            if ~isempty(prod.replacefds1)
                tbl1 = replacefield(tbl1, prod.replacefds1);
            end
            if ~isempty(prod.replacefds2)
                tbl2 = replacefield(tbl2, prod.replacefds2);
            end

            % some sanity checks on the table's field names.
            fds1 = tbl1.fdnames;
            fds2 = tbl2.fdnames;

            assert(all(ismember(mergefds, fds1)), ['There exist merge fields that do ' ...
                                'not belong to fields of first table']);
            assert(all(ismember(mergefds, fds2)), ['There exist merge fields that do ' ...
                                'not belong to fields of second table']);
            ofds1 = fds1(~ismember(fds1, crossfds));
            ofds2 = fds2(~ismember(fds2, crossfds));
            assert(all(~ismember(ofds1, ofds2)) & all(~ismember(ofds1, ofds2)), ...
                   ['There exist fields with same name in first and second ' ...
                    'table that are neither merged or reduced.']);
            
            if isempty(prod.pivottbl)
                
                [pivottbl, indstruct] = crossIndexArray(tbl1, tbl2, crossfds);
                
                dispind1 = indstruct{1}.inds;
                dispind2 = indstruct{2}.inds;
                
            else

                pivottbl = prod.pivottbl;

                pivotfds = pivottbl.fdnames;
                
                mergefds1 = tbl1.fdnames(ismember(tbl1.fdnames, pivotfds));

                map = TensorMap();
                map.fromTbl  = tbl1;
                map.toTbl    = pivottbl;
                map.mergefds = mergefds1;
                dispind1 = getDispatchInd(map);
                
                mergefds2 = tbl2.fdnames(ismember(tbl2.fdnames, pivotfds));

                map = TensorMap();
                map.fromTbl  = tbl2;
                map.toTbl    = pivottbl;
                map.mergefds = mergefds2;
                dispind2 = getDispatchInd(map);
                
            end
            
            if isempty(prod.tbl3)
                fds1 = tbl1.fdnames;
                fds1 = fds1(~ismember(fds1, horzcat(crossfds, reducefds1)));
                fds2 = tbl2.fdnames;
                fds2 = fds2(~ismember(fds2, horzcat(crossfds, reducefds2)));
                fds3 = {fds1{:}, fds2{:}, mergefds{:}};
                [tbl3, dispind3] = projIndexArray(pivottbl, fds3);
                dispind3 = dispind3.inds;
            else
                
                tbl3 = prod.tbl3;
                fds3 = tbl3.fdnames;
                
                fds3 = fds3(~ismember(fds3, reducefds1));
                fds3 = fds3(~ismember(fds3, reducefds2));

                if prod.needIntersectTbl3

                    [pivottbl, indstruct] = crossIndexArray(pivottbl, tbl3, fds3);
                    disppivot = indstruct{1}.inds;

                    dispind1 = dispind1(disppivot);
                    dispind2 = dispind2(disppivot);
                    
                end
                
                map = TensorMap();
                map.fromTbl = tbl3;
                map.toTbl = pivottbl;
                map.mergefds = fds3;
                dispind3 = getDispatchInd(map);
                
            end
            
            if prod.settbl3
                replacefds3 = prod.replacefds3;
                if ~isempty(replacefds3)
                    tbl3 = replacefield(tbl3, replacefds3);
                end
                prod.tbl3 = tbl3;
            end

            prod.pivottbl   = pivottbl;
            
            prod.dispind1  = dispind1;
            prod.dispind2  = dispind2;
            prod.dispind3  = dispind3;
            
            prod.issetup = true;
            
        end

        function prodAB = eval(prod, A, B)
            assert(prod.issetup, ['tensor product is not setup. Use method ' ...
                                'setup']);

            if prod.useMex

                prodAB = tensorProduct(prod.tbl3.num        , ...
                                       double(A)            , ...
                                       double(B)            , ...
                                       uint64(prod.dispind1), ...
                                       uint64(prod.dispind2), ...
                                       uint64(prod.dispind3));
                return
                
            else
                
                chunksize = prod.chunksize;
                
                dispind1 = prod.dispind1;
                dispind2 = prod.dispind2;
                dispind3 = prod.dispind3;
                
                n3 = prod.tbl3.num;
                n  = prod.pivottbl.num;

                if ~isempty(chunksize) && isa(A, 'double') && isa(B, 'double')

                    pivotsize = numel(dispind1); % size is same for dispind1, dispind2 and dispind3
                    nchunks   = ceil(pivotsize/chunksize);

                    if mrstVerbose() > 0
                        fprintf('number of chunks %d ', nchunks);
                    end
                    
                    prodAB = zeros(n3, 1, 'double');
                    
                    for ichunk = 1 : nchunks

                        if ichunk < nchunks
                            ind = (1 + (ichunk - 1)*chunksize) : ichunk*chunksize;
                            if ichunk == 1
                                prodABc = zeros(chunksize, 1);
                                pind = (1 : chunksize)';
                            end
                        else
                            npind = pivotsize - (nchunks - 1)*chunksize;
                            prodABc = zeros(npind, 1);
                            pind = (1 : npind)';
                            ind = (1 + (ichunk - 1)*chunksize) : pivotsize;
                        end

                        prodABc(pind) = A(dispind1(ind)).*B(dispind2(ind));
                        
                        prodAB = prodAB + accumarray(dispind3(ind), prodABc, [n3, 1]);

                        if mrstVerbose() > 0
                            fprintf('.');
                        end
                        
                    end
                    
                    if mrstVerbose() > 0
                        fprintf('\n');
                    end
                    
                else
                    
                    if isa(A, 'double') && isa(B, 'double')
                        
                        prodAB = accumarray(dispind3, A(dispind1).*B(dispind2), [n3, 1]);
                        
                    else

                        prodAB = A(dispind1).*B(dispind2);
                        M = sparse(dispind3, (1 : n)', 1, n3, n);
                        prodAB = M*prodAB;
                        
                    end
                    
                end

            end
            
        end

        function isok = checkSetup(prod, varargin)
        % check if the dispatching indices (dispinds) are setup correctly with respect to the pivot IndexArray (prod.pivottbl)
            
            pivottbl = prod.pivottbl;
            tbl1 = prod.tbl1;
            tbl2 = prod.tbl2;
            tbl3 = prod.tbl3;
            
            if ~isempty(prod.replacefds1)
                tbl1 = replacefield(tbl1, prod.replacefds1);
            end
            
            if ~isempty(prod.replacefds2)
                tbl2 = replacefield(tbl2, prod.replacefds2);
            end
            
            if ~isempty(prod.replacefds3)
                tbl3 = replacefield(tbl3, prod.replacefds3);
            end
            
            if (nargin > 0) & ~isempty(varargin)
                pivottbl = replacefield(pivottbl, varargin{1});
            end
            
            map = TensorMap();
            map.fromTbl = tbl1;
            map.toTbl = pivottbl;
            map.mergefds = tbl1.fdnames;
            dispind1 = map.getDispatchInd();
            
            isok = all(dispind1 == prod.dispind1);
            
            map = TensorMap();
            map.fromTbl = tbl2;
            map.toTbl = pivottbl;
            map.mergefds = tbl2.fdnames;
            dispind2 = map.getDispatchInd();
            
            isok = isok & all(dispind2 == prod.dispind2);
            
            map = TensorMap();
            map.fromTbl = tbl3;
            map.toTbl = pivottbl;
            map.mergefds = tbl3.fdnames;
            dispind3 = map.getDispatchInd();
            
            isok = isok & all(dispind3 == prod.dispind3);
            
            
        end

        function M = setupMatrix(prod, vals)

            M = SparseTensor();
            M = M.setFromTensorProd(vals, prod);
            M = M.getMatrix();

        end

        
        function matrices = getMatrices(prod)
            
            dispind1 = prod.dispind1;
            dispind2 = prod.dispind2;
            dispind3 = prod.dispind3;
            
            n1 = prod.tbl1.num;
            n2 = prod.tbl2.num;
            n3 = prod.tbl3.num;

            npiv = prod.pivottbl.num;

            % Dispatch from tbl1 to pivot
            D1 = sparse((1 : npiv)', dispind1, 1, npiv, n1);
            % Dispatch from tbl2 to pivot
            D2 = sparse((1 : npiv)', dispind2, 1, npiv, n2);
            % Sum from pivot to tbl3
            S = sparse(dispind3, (1 : npiv)', 1, n3, npiv);

            matrices = struct('D1', D1, ...
                              'D2', D2, ...
                              'S', S);
            
        end
            
        function [ind1, ind2] = getDispatchInd(prod)
        % In the case where the product is set up to create a bilinear mapping (see
        % SparseTensor.setFromTensorProd method), then we can use this function
        % to obtain the dispatching indices. It can be only used if the pivot
        % space (given by pivottbl) is the same as the space given by
        % tbl1. Hence, the assert statement below, which should be enough to
        % detect those cases.
            
            tbl1 = prod.tbl1;
            tbl2 = prod.tbl2;
            tbl3 = prod.tbl3;
            
            if ~isempty(prod.replacefds1)
                tbl1 = replacefield(tbl1, prod.replacefds1);
            end
            if ~isempty(prod.replacefds2)
                tbl2 = replacefield(tbl2, prod.replacefds2);
            end
            
            isok = all(ismember(tbl2.fdnames, tbl1.fdnames));
            isok = isok & all(ismember(tbl3.fdnames, tbl1.fdnames));
            assert(isok, ['getDispatchInd cannot be used for this table (not ' ...
                          'well defined in this case)']);
            
            map = TensorMap();
            map.fromTbl = tbl3;
            map.toTbl = tbl1;
            map.mergefds = tbl3.fdnames;
            ind1 = getDispatchInd(map);
            
            map = TensorMap();
            map.fromTbl = tbl2;
            map.toTbl = tbl1;
            map.mergefds = tbl2.fdnames;
            ind2 = getDispatchInd(map);
            
        end
        
        
    end
    
        
   
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
