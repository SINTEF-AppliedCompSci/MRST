classdef TensorMap 
%
%
% SYNOPSIS:
%   TensorMap(varargin)
%
% DESCRIPTION:
% Given two IndexArrays (fromTbl and toTbl) and an indexing described by its
% name (mergefds), there exists a unique mapping from fromTbl to toTbl. This
% class sets up this mapping. We use is as follows: First, we create an instance
% (map = TensorMap), then we assign the argument and result space for the
% mapping (assign fromTbl, toTbl) and the field names that will be merged by
% intersection, see explanation below (we assign mergefds). Finally, we setup
% the mapping by running map = map.setup(). In this last stage, an internal
% representation of the mapping (using local linear indices) is created and the
% mapping of vectors using (v = map.eval(u)) will be computed from that
% efficiently. This representation is given by the index vectors dispind1
% dispind2 and the IndexArray pivottbl. In particular, the sparsity for an
% efficient computation of the product is computed (see paper on Tensor Product
% assembly).
%
% Let us explain how the mapping (denoted map) is defined. We consider
%
%  * u[i, j] where the index [i, j] belongs to fromTbl
%  * v[i, k] where the index [i, j] belongs to toTbl
%
% where i corresponds to the indexing mergedfs. Note that i, j, k can be
% themselves multiple-indices.
%
% We have v = map(u) if
%
%   v[i, k] = sum_{j} u[i, j]
%
% We consider only i such that [i, k] belongs to toTbl and, in the summation we
% consider only j such that [i, j] belongs to fromTbl. To do so, we have thus
% computed an kind of intersection between fromTbl and toTbl, taking all the
% index i of the indexing mergefds that appear in both fromTbl and toTbl
% 
% In the setup of the mapping (before calling map = map.setup()), the fieldnames
% of fromTbl and toTbl can be changed using prod.replaceFromTblfds and
% prod.replaceToTblfds. This option is used very often as different names for
% the same indexing appear naturally.
%
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
%   `IndexArray`.
    
    properties
        fromTbl           % IndexArray for the argument space of the mapping
        toTbl             % IndexArray for the output space of the mapping
        mergefds          % Field names for merging
         
        replaceFromTblfds % Possibility to change field names of first table
                          % (before setting up map)
        replaceToTblfds   % Possibility to change field names of second table
                          % (before setting up map)
        
        % Structures used for the computation. Initialized only once.
        pivottbl          % Structure for the pivot tbl which is used to setup the
                          % map (note that this structure is not necessary to be set
                          % for the map to be setup. It is enough to have the
                          % dispatch indices below)
        
        % dispatch indices 
        dispind1          % from fromTbl to pivottbl 
        dispind2          % from toTbl to pivottbl
        
        issetup           % Flag is set to true is map has been set up.
    end
   
    methods
        
        function map = TensorMap(varargin)
            
            map.fromTbl = [];
            map.toTbl = [];
            map.replaceFromTblfds = [];
            map.replaceToTblfds = [];
            map.mergefds = [];
            map.pivottbl = [];
            map.dispind1 = [];
            map.dispind2 = [];
            map.issetup = false;
            
            map = merge_options(map, varargin{:}); 
            
            if isempty(map.mergefds)
                map.mergefds = {};
            end            
            
        end
        
        function map = setup(map)
        % Compute dispind1, dispind2 and pivottbl
            
            if map.issetup
                fprintf(['TensorMap instance already setup. I am not changing', ...
                         ' it.\n']);
                return
            end
            
            mergefds  = map.mergefds;
            fromTbl = map.fromTbl;
            toTbl = map.toTbl;
            
            if ~isempty(map.replaceFromTblfds)
                fromTbl = replacefield(fromTbl, map.replaceFromTblfds);
            end
            if ~isempty(map.replaceToTblfds)
                toTbl = replacefield(toTbl, map.replaceToTblfds);
            end

            % some sanity checks on the table's field names.
            fds1 = fromTbl.fdnames;
            fds2 = toTbl.fdnames;
            assert(numel(unique(fds1)) == numel(fds1), ['repetition in first ' ...
                                'table field names']);
            assert(numel(unique(fds2)) == numel(fds2), ['repetition in second ' ...
                                'table field names']);            
            assert(all(ismember(mergefds, fds1)), ['merge  fields do ' ...
                                'not belong to fields of first table']);
            assert(all(ismember(mergefds, fds2)), ['merge  fields do ' ...
                                'not belong to fields of second table']);
            assert(all(ismember(mergefds, fds1) | ismember(mergefds, fds2)), ['replace  fields do ' ...
                                'not belong to fields of either first or second table']);
            ofds1 = fds1(~ismember(fds1, mergefds));
            ofds2 = fds2(~ismember(fds2, mergefds));
            assert(all(~ismember(ofds1, ofds2)) & all(~ismember(ofds1, ofds2)), ...
                   ['There exist fields with same name in first and second ' ...
                    'table that are neither merged or reduced.']);
            
            if isempty(map.pivottbl)
                [pivottbl, indstruct] = crossIndexArray(fromTbl, toTbl, mergefds);
                dispind1 = indstruct{1}.inds;
                dispind2 = indstruct{2}.inds;
                map.pivottbl = pivottbl;
            else
                givenpivottbl = map.pivottbl;
                [pivottbl, indstruct] = crossIndexArray(fromTbl, toTbl, mergefds);                
                dispind1 = instruct{1}.inds;
                dispind2 = instruct{2}.inds;
                % pivottbl and givenpivottbl corresponds to the same multi-indices
                % but not necessarily in the same order. 
                % We check that they at least have the same field
                givenfds = givenpivottbl.fdnames;
                fds = pivottbl.fdnames;
                assert(all(ismember(fds, givenfds)) & all(ismember(givenfds, ...
                                                                  fds)), ...
                       'non matching fields');
                [~, indstruct] = crossIndexArray(pivotbl, givenpivotbl, givenfds);
                ind1 = indstruct{1}.inds;
                ind2 = indstruct{2}.inds;
                indi2(ind2) = (1 : givenpivottbl.num)';
                ind = ind1(ind2);
                dispind1 = dispind1(ind);
                dispind2 = dispind2(ind);
            end
            
            map.dispind1  = dispind1;
            map.dispind2  = dispind2;
            
            map.issetup = true;
            
        end

        function v = eval(map, u)
        % Evaluate mapping to a vector u, which follows indexing of fromTbl
            assert(map.issetup, ['tensor map is not setup. Use method ' ...
                                'setup']);
            
            dispind1 = map.dispind1;
            dispind2 = map.dispind2;
            toTbl    = map.toTbl;
            fromTbl  = map.fromTbl;
            
            if isa(u, 'ADI')
                M = sparse(dispind2, dispind1, 1, toTbl.num, fromTbl.num);
                v = M*u;
            else
                v = u(dispind1);
                v = accumarray(dispind2, v, [toTbl.num, 1]);
            end
            
        end
        
        function inds = getDispatchInd(map)
        % Compute dispatching indices: For an index i in map.toTbl gives the index number in map.fromTbl that will be sent to i
        % by the mapping map
            
            map = map.setup();
            
            ind1 = map.dispind1;
            ind2 = map.dispind2;

            toTbl = map.toTbl;
            
            indi2 = zeros(toTbl.num, 1);
            indi2(ind2) = (1 : toTbl.num)';
            inds = ind1(indi2);
            
        end
        
        function v = evalTranspose(map, u)
        % Evaluate the transpose of the mapping
            assert(map.issetup, ['tensor map is not setup. Use method ' ...
                                'setup']);
            
            dispind1 = map.dispind1;
            dispind2 = map.dispind2;
            
            fromTbl = map.fromTbl;
            
            v = u(dispind2);
            v = accumarray(dispind1, v, [fromTbl.num, 1]);
            
        end

        function isok = checkSetup(map, varargin)
        % check if the dispatching indices (dispinds) are setup correctly with respect to the pivot IndexArray (prod.pivottbl)
            
            pivottbl = map.pivottbl;
            fromtbl  = map.fromTbl;
            totbl    = map.toTbl;
            
            if ~isempty(map.replaceFromTblfds)
                fromtbl = replacefield(fromtbl, map.replaceFromTblfds);
            end
            
            if ~isempty(map.replaceToTblfds)
                totbl = replacefield(totbl, map.replaceToTblfds);
            end
            
            if (nargin > 0) & ~isempty(varargin)
                pivottbl = replacefield(pivottbl, varargin{1});
            end
            
            dismap = TensorMap();
            dismap.fromTbl  = fromtbl;
            dismap.toTbl    = pivottbl;
            dismap.mergefds = fromtbl.fdnames;
            dispind1 = dismap.getDispatchInd();
            
            isok = all(dispind1 == map.dispind1);
            
            dismap = TensorMap();
            dismap.fromTbl  = totbl;
            dismap.toTbl    = pivottbl;
            dismap.mergefds = totbl.fdnames;
            dispind2 = dismap.getDispatchInd();
            
            isok = isok & all(dispind2 == map.dispind2);
            
        end

        function mat = getMatrix(map)
        % get matlab sparse matrix representation of the mapping in term of the local
        % indices of fromTbl and toTbl.
            n1 = map.fromTbl.num;
            n2 = map.toTbl.num;
            n  = map.pivottbl.num;
            
            ind1 = map.dispind1;
            ind2 = map.dispind2;
            
            mat1 = sparse((1 : n)', ind1, 1, n, n1);
            mat2 = sparse((1 : n)', ind2, 1, n, n2);
            
            mat = mat2'*mat1;
            
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
