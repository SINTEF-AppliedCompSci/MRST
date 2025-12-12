classdef SparseTensorMap

    
    properties
        
        tensor
        
        toTbl             % IndexArray for the output space of the mapping
        mergefds          % Field names for merging
         
        replaceFromTblfds % Possibility to change field names of first table
                          % (before setting up map)
        replaceToTblfds   % Possibility to change field names of second table
                          % (before setting up map)
        
    end
   
    methods
        
        function map = SparseTensorMap(tensor, varargin)

            opts = struct('toTbl'            , [], ...
                          'mergefds'         , [], ...
                          'replaceFromTblfds', [], ...
                          'replaceToTblfds'  , []);
            
            map = merge_options(map, varargin{:}); 
            
            if isempty(map.replaceFromTblfds)
                map.replaceFromTblfds = {};
            end            

            if isempty(map.replaceToTblfds)
                map.replaceToTblfds = {};
            end            

            if isempty(map.mergefds)
                map.mergefds = {};
            end            
            
        end

        function B = eval(map)

            mergefds = map.mergefds;
            fromTbl  = map.tensor.tbl;
            toTbl    = map.toTbl;
            
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

            lA = fromTbl.gets(mergefds);
            jA = fromTbl.gets(reducefds);

            lB = toTbl.gets(mergefds);
            kB = toTbl.gets(ofds2);

            A = map.tensor.vals;

            [m_l, ~, b_l] = unique([lA; lB], 'rows');
            b_lA = b_l(1 : size(lA, 1));
            b_lB = b_l(size(lA, 1) + (1 : size(lB, 1)));

            [m_i, ~, b_iA] = unique(iA, 'rows');
            [m_j, ~, b_jA] = unique(jA, 'rows');
            [m_k, ~, b_kB] = unique(kB, 'rows');

            b_l = intersect(b_lA, b_lB, 'sorted');
                        
            LIA = ismember(b_lA, b_l);
            b_lA = b_lA(LIA);
            b_iA = b_iA(LIA);
            b_jA = b_jA(LIA);
            A    = A(LIA);

            LIA = ismember(b_lB, b_l);
            b_lB = b_lB(LIA);
            b_kB = b_kB(LIA);

            [b_lA, isort] = sort(b_lA);
            b_iA = b_iA(isort);
            b_jA = b_jA(isort);
            A = A(isort);
 
            [b_lB, isort] = sort(b_lB);
            b_kB = b_kB(isort);

            b_lA_b_iA = unique([b_lA, b_iA], 'rows');
            [~, nI] = rlencode(b_lA_b_iA(: , 1));

            b_lB_b_kB = unique([b_lB, b_kB], 'rows');
            [~, nK] = rlencode(b_lB_b_kB(: , 1));

            nalloc = sum(nI.*nK);

            [u_b_lA, n_lA] = rlencode(b_lA);
            [u_b_lB, n_lB] = rlencode(b_lB);
            
            assert(all(u_b_lA == u_b_lB));
            u_b_l = u_b_lA;
            
            startA   = 1;
            endA     = n_lA(1);
            
            startB   = 1;
            endB     = n_lB(1);
            
            b_iC_all = zeros(nalloc, 1, 'uint64');
            b_kC_all = zeros(nalloc, 1, 'uint64');
            B_all    = zeros(nalloc, 1, 'double');
            
            n_u_b_l = [];
            
            posB = 0;

            for ind_u_b_l = 1 : numel(u_b_l)
                
                ind = (startA : endA);
                sA = sparse(b_iA(ind), b_jA(ind), A(ind));
                sB = sum(sA, 2)
                
                [b_iB, ~, B] = find(sB);                    

                nB = size(B, 1);
                
                n_u_b_l(end + 1) = nB;
                
                b_iB_all(posB  + (1 : nB)) =  b_iB;
                B_all(posB  + (1 : nB))    =  B ;
                
                posB = posB + nB;

                if ind_u_b_l < numel(u_b_l)
                    startA = endA + 1;
                    endA   = startA + n_lA(ind_u_b_l + 1) - 1;                        
                end

            end
    
            fprintf('nalloc : %d, posC : %d, difference : %d\n', nalloc, posB, abs(nalloc - posB));
            
            B = B_all(1 : posB);

            lB_all = rldecode(m_l(u_b_l, :), n_u_b_l);
            fdnames = {mergefds{:}, ofds1{:}, ofds2{:}};
            tbl = IndexArray([], 'fdnames', fdnames, 'inds', [lB_all, m_i(b_iB_all(1 : posB), :), m_k(b_kB_all(1 : posB), :)]);

            B = SparseTensor(B, tbl);

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
