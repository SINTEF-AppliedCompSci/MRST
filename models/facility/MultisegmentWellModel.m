classdef MultisegmentWellModel < WrapperModel
    
    properties
        radius = 9*centi*meter
        trajectory
    end
    
    methods
   
        function model = MultisegmentWellModel(rmodel, coords, varargin)
            
            rmodel = removeStateFunctionGroupings(rmodel);
            model = model@WrapperModel(rmodel);
            trajectory = computeTraversedCells(rmodel.G, coords);
            trajectory.coords = coords;
            
            model.trajectory = trajectory;
            
            model = merge_options(model, varargin{:});
            
            G = extractSubgrid(rmodel.G, trajectory.cell);
            rock = extractSubRock(rmodel.G, rmodel.rock, trajectory.cell);
            
            model.G = G;
            model.parentModel.G = G;
            model.parentModel.rock = rock;
            model = model.setupOperators(G, rock);
            
        end
        
        function model = setupOperators(model, G, rock, varargin)
            % Well segment volumes
            model.parentModel = model.parentModel.setupOperators(G, rock);
            length = sqrt(sum(model.trajectory.vec.^2, 2));
            area   = pi*model.radius.^2;
            vol    = length.*area;
            model.operators.pv = vol;
            % Well-segment-to-matrix transmissibility
            model.parentModel.operators.WI = computeWellIndex(G, rock, model.radius, (1:G.cells.num)');
            % TODO: Well-segmet-to-well-segment transmissibility
            
        end
        
        function [p, state] = getProp(model, state, name)
            [p, state] = model.parentModel.getProp(state, name);
        end
         
        function [eqs, names, types, state] = getModelEquations(model, state0, state, dt, drivingForces)
             [eqs, names, types, state] = model.parentModel.getModelEquations(state0, state, dt, drivingForces);
         end
        
    end
    
end

function rock = extractSubRock(G, rock, cells)

    names = fieldnames(rock)';
    
    for name = names
        v = rock.(name{1});
        if size(v,1) ~= G.cells.num, continue; end
        rock.(name{1}) = v(cells,:);
    end

end


%-------------------------------------------------------------------------%
function cells = topoSortWellCells(G, cells, reverse)
% Topologically sort well cells. Mostly for visualization - should not have
% any effect on drawdown computations.

    n = numel(cells);
    if n == 1, return; end
    
    if nargin < 3, reverse = false; end
    cells0 = cells; %#ok For debugging
    
    loc2glob = cells;
    glob2loc = nan(G.cells.num, 1); glob2loc(cells) = (1:n)';

    facePos = mcolon(G.cells.facePos(cells), G.cells.facePos(cells+1)-1);
    faces   = G.cells.faces(facePos,1);
    N       = G.faces.neighbors(faces,:);
    keep = all(ismember(N, cells), 2);
    N = N(keep,:); 
    N = glob2loc(N);
    N = unique(N, 'rows');
    N = [N; fliplr(N)];
    A = sparse(N(:,1), N(:,2), 1, n, n);
   
    endpoints = find(sum(A,2) == 1);
    
    [~, ix] = min(G.cells.centroids(endpoints,3));
    first = endpoints(ix);
 
    cells = nan(n,1);
    current = false(n,1);
    current(first) = true;
    
    processed = false(n,1);
    pno = 1;
    while any(~processed)
        cells(pno) = find(current);
        processed(current) = true;
        candidates = A*current > 0 & ~processed;
        if nnz(candidates) > 1
            warning('Branching wells not supported')
            cix = find(candidates);
            candidates(cix(2:end)) = false;
            processed(cix(2:end)) = true;
        end
        
        current = candidates;
        pno = pno + 1;
    end
    cells = cells(~isnan(cells));
    cells = loc2glob(cells);
    
    if reverse, cells = flip(cells); end

end