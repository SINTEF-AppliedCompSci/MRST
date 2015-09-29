run('../../matlab/project-mechanics-fractures/mystartup.m')

                            %   Build cartesian grid and compute geometry.
nx = 3; ny = 3;             %   Grid dimensions.
G = cartGrid([nx, ny]);
G = mrstGridWithFullMappings(G);
G = computeGeometry(G);
n = G.cells.num;            %   Number of cells in G.

pos = 0;
for c = 1:n
    
    nodeNum = G.cells.nodePos(c):G.cells.nodePos(c+1)-1
    nodes = G.cells.nodes(nodeNum)
    X = G.nodes.coords(nodes,:)
    Sl = locS(X);
    
    dofvec=mcolon(2*(nodes-1)+1,2*(nodes-1)+2);
    iindl=repmat(dofvec',1,numel(dofvec));
    jindl=repmat(dofvec,numel(dofvec),1);
    nnzl = numel(iindl);
    iglob(pos+1:pos+nnzl) = iindl;
    jglob(pos+1:pos+nnzl) = jindl;
    mglob(pos+1:pos+nnzl) = Sl;
    pos = pos+nnzl;
end
S=sparse(iglob(:),jglob(:),mglob(:),ndof,ndof)