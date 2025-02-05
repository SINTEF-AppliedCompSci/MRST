classdef CrossIndexArrayGenerator 

    properties
        
        tbl1        % IndexArray for the first table
        tbl2        % IndexArray for the second table
        mergefds    % Field names for merging

        replacefds1 % Possibility to change field names of first table
                    % (before setting up gen)
        replacefds2 % Possibility to change field names of second table
                    % (before setting up gen)

        sortfds                   % Field names for sorting (if empty, no sorting is done)
        sortKeepAllFields = true; % If true, all fields are kept in the sorted table
        
        %% dispactch mappings updated during setup
        tbl
        ind1
        ind2

        %% options
        opts = {};
        
    end

    methods

        function gen = CrossIndexArrayGenerator(varargin)

            gen.tbl1        = [];
            gen.tbl2        = [];
            gen.replacefds1 = [];
            gen.replacefds2 = [];
            gen.mergefds    = [];
            gen.sortfds     = [];
            
            gen = merge_options(gen, varargin{:}); 

            if isempty(gen.mergefds)
                gen.mergefds = {};
            end            

        end

        function [tbl, gen] = eval(gen)


            mergefds  = gen.mergefds;
            tbl1 = gen.tbl1;
            tbl2 = gen.tbl2;

            if ~isempty(gen.replacefds1)
                 tbl1 = replacefield(tbl1, gen.replacefds1);
            end

            if ~isempty(gen.replacefds2)
                 tbl2 = replacefield(tbl2, gen.replacefds2);
            end

            [tbl, indstruct] = crossIndexArray(tbl1, tbl2, mergefds, gen.opts{:});

            gen.tbl = tbl;
            gen.ind1 = indstruct{1}.inds;
            gen.ind2 = indstruct{2}.inds;
            
            if ~isempty(gen.sortfds)
                [tbl, dispind] = sortIndexArray(tbl, gen.sortfds, 'keepAllFields', gen.sortKeepAllFields);
                gen.tbl = tbl;
                gen.ind1 = gen.ind1(dispind);
                gen.ind2 = gen.ind2(dispind);
            end
            
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
