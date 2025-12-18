classdef SparseTensorMap

    
    properties (SetAccess = immutable)

        tensor

    end

    properties
        
        toTbl             % IndexArray for the output space of the mapping
        mergefds          % Field names for merging
         
        replaceFromTblfds % Possibility to change field names of first table
                          % (before setting up map)
        replaceToTblfds   % Possibility to change field names of second table
                          % (before setting up map)
        
    end
   
    methods
        
        function map = SparseTensorMap(tensor, varargin)

            map = merge_options(map, varargin{:}); 

            map.tensor = tensor;
            
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
            reducefds = fds1(~ismember(fds1, mergefds));
            ofds = fds2(~ismember(fds2, mergefds));
            assert(all(~ismember(reducefds, ofds)) & all(~ismember(reducefds, ofds)), ...
                   ['There exist fields with same name in first and second ' ...
                    'table that are neither merged or reduced.']);

            lA = fromTbl.gets(mergefds);
            iA = fromTbl.gets(reducefds);

            lB = toTbl.gets(mergefds);

            A = map.tensor.vals;

            [m_l, ~, b_l] = unique([lA; lB], 'rows');
            b_lA = b_l(1 : size(lA, 1));
            b_lB = b_l(size(lA, 1) + (1 : size(lB, 1)));

            [m_i, ~, b_iA] = unique(iA, 'rows');

            b_l = intersect(b_lA, b_lB, 'sorted');
                        
            LIA = ismember(b_lA, b_l);
            b_lA = b_lA(LIA);
            b_iA = b_iA(LIA);
            A    = A(LIA);

            sA = sparse(b_lA, b_iA, A);
            sB = sum(sA, 2);

            B = zeros(toTbl.num, 1);

            [~, LOCB] = ismember(b_lB, b_lA);

            indB  = find(LOCB);
            indsB = b_lA(LOCB(indB));
            
            B(indB) = sB(indsB);

            B = SparseTensor(B, toTbl);

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
