function [I_basis, I, bnd] = miteratedJacobiBasis(A, CG, varargin)
    % Force A to be M-matrix-like
    A = A - diag(sum(A, 2));
    
    Ad = diag(A);
    [ii, jj, vv] = find(A - diag(diag(A)));
    
    assert(~any(ii == jj));
    assert(issorted(jj));
    
    [a, b] = rlencode(jj);
    jpos = [0; cumsum(b)] + 1;
    
    G = CG.parent;

    % Hacking to get the boundary regions via some intersections etc
    bnd = cell(CG.cells.num, 1);
    
    interactionCount = ones(G.cells.num, 1);
    for i = 1:CG.cells.num
        interactionCount(CG.cells.interaction{i}) =...
                        interactionCount(CG.cells.interaction{i}) + 1;
    end
    notinter = find(interactionCount > 1);
    
    for i = 1:CG.cells.num
        c = bndregion(G, CG.cells.interaction{i});
        c = intersect(c, notinter);
        bnd{i} = c;
    end
    
    % The boundary of each interaction region
    bndreg = vertcat(bnd{:});
    bndregPos = [0; cumsum(cellfun(@numel, bnd))] + 1;
    
    otherEdges = cell(CG.cells.num, 1);
    for i = 1:CG.cells.num
        otherEdges{i} = intersect(CG.cells.interaction{i}, bndreg);
    end
    otheredge = vertcat(otherEdges{:});
    otheredgePos = [0; cumsum(cellfun(@numel, otherEdges))] + 1;
    
    % Interaction region for each coarse block
    interact = vertcat(CG.cells.interaction{:});
    numint = cellfun(@numel, CG.cells.interaction);
    interactPos = [0; cumsum(numint)] + 1;
    interactCoarse = rldecode((1:CG.cells.num)', numint);
    
    % Initial guess
    I0 = 0*interact;
    for i = 1:CG.cells.num
        pos = interactPos(i):(interactPos(i+1)-1);
        I0(pos) = ismember(interact(pos), find(CG.partition == i));
    end
    
    % The overlap - regions where normalization is required
    overlap = unique(interact);
    
    ii = int32(ii) - 1;
    jpos = int32(jpos) - 1;
    bndreg = int32(bndreg) - 1;
    bndregPos = int32(bndregPos) - 1;
    interact = int32(interact) - 1;
    interactPos = int32(interactPos) - 1;
    overlap = int32(overlap) - 1;
    
    otheredgePos = int32(otheredgePos) - 1;
    otheredge = int32(otheredge) - 1;
    
    disp('Starting compute in compiled code!');
    tic()
    I = mex_iteratedJacobiBasis(I0, int32(numel(I0)),...
                            int32(G.cells.num), int32(CG.cells.num), ...
                            Ad, ii, jpos, vv, ...
                            bndreg, bndregPos, ...
                            interact, interactPos, ...
                            otheredge, otheredgePos, ...
                            overlap, int32(numel(overlap)));
    toc()
%                                 neigh, neighPos, ...

%     I = I0;
    
    I_basis = sparse(double(interact) + 1, interactCoarse, I, G.cells.num, CG.cells.num);
    
end


function c = bndregion(g, sub)
    if islogical(sub)
        sub = find(sub);
    end


   present          = false([g.cells.num + 1, 1]);
   present(sub + 1) = true;
   present( 0  + 1) = true;

   tmp = present(g.faces.neighbors + 1);
   f = xor(tmp(:,1), tmp(:,2));


  n = double(g.faces.neighbors(f,:));
  c = sum(n .* double(present(n + 1)), 2);

  c = c(c~=0);
end
