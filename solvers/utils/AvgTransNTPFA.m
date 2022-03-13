function T = AvgTransNTPFA(G, OSflux)
%Undocumented Utility Function

%{
Copyright 2009-2022 SINTEF Digital, Mathematics & Cybernetics.

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

    dispif(mrstVerbose, 'AvgTransNTPFA... ');
    timer = tic;

    nf = G.faces.num;
    nc = G.cells.num;

    celltbl.cells = (1 : nc)';
    celltbl = IndexArray(celltbl);
    
    facetbl.faces = (1 : nf)';
    facetbl = IndexArray(facetbl);
    
    internal = (1 : nf)';
    internal(~all(G.faces.neighbors ~= 0, 2)) = [];
    intfacetbl.faces = internal;
    intfacetbl = IndexArray(intfacetbl);
    
    T = cell(2, 1);
    
    % We unpack the OSflux to multi-index format using IndexArray
    cells   = cell(2, 1);
    faces   = cell(2, 1);
    weights = cell(2, 1);
        
    for j = 1 : 2
        osflux = OSflux(:, j);
        % We remove the last entry in each cell (equal to zero, I do not know why.)
        osflux = cellfun(@(x) x(1 : (end - 1), :), osflux, 'uniformoutput', false);
        ossizes = cellfun(@(x) size(x), osflux, 'uniformoutput', false);
        osflux = cell2mat(osflux);
        cells{j} = osflux(:, 1);
        weights{j} = osflux(:, 2);
            
        ossizes = cell2mat(ossizes);
        faces{j} = rldecode((1 : nf)', ossizes(:, 1));
    end
    
    lincellfacetbls = cell(2, 1);
    
    for j = 1 : 2
        
        clear lincellfacetbl
        lincellfacetbl.cells = cells{j};
        lincellfacetbl.faces = faces{j};
        lincellfacetbl = IndexArray(lincellfacetbl);
        
        % We only consider the internal faces
        linintcellfacetbl = crossIndexArray(lincellfacetbl, intfacetbl , {'faces'});
        
        map = TensorMap();
        map.fromTbl = lincellfacetbl;
        map.toTbl = linintcellfacetbl;
        map.mergefds = {'cells', 'faces'};
        map = map.setup();
        
        weights{j} = map.eval(weights{j});
        lincellfacetbls{j} = linintcellfacetbl;
        
    end
    
    %% We have to change the weights of the "inward" cell to opposite sign
    
    intfaces = intfacetbl.get('faces');
    
    for j = 1 : 2
        
        clear cellfacetbl;
        cellfacetbl.faces = intfaces;
        cellfacetbl.cells = G.faces.neighbors(intfaces, j);
        cellfacetbl = IndexArray(cellfacetbl);

        map = TensorMap();
        map.fromTbl = cellfacetbl;
        map.toTbl = lincellfacetbls{j};
        map.mergefds = {'cells', 'faces'};
        map = map.setup();
        
        chsign = logical(map.eval(ones(cellfacetbl.num, 1)));
        
        weights{j}(chsign) = -weights{j}(chsign);
        
    end        
    
    %% We setup the flux mappings
    
    for j = 1 : 2
        
        prod = TensorProd();
        prod.tbl1 = lincellfacetbls{j};
        prod.tbl2 = celltbl;
        prod.tbl3 = intfacetbl;
        prod.reducefds = {'cells'};
        prod = prod.setup();
        
        T{j} = SparseTensor();
        T{j} = T{j}.setFromTensorProd(weights{j}, prod);
        
        T{j} = T{j}.getMatrix();
        
    end
    
end

