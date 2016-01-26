function I = polyhedronInt(G,f)

%   http://people.sc.fsu.edu/~jburkardt/m_src/tetrahedron_arbq_rule/tetrahedron_arbq_rule.html

Xq = [-.2677028585910073,  .5939017006199488E-01, -.3969426941142150; ...
       .1510931841533142,  .4354087732476309,      .2151719254306919; ...
      -.1367699195237390, -.3280341115590410,      .2950846519937133; ...
       .8067449309051964, -.3400977314285288,     -.3430378951002550];

w  = [ 0.2918008865477151, 0.2706884392356724, ...
      0.3098349601753470, 0.9865925745591227E-01 ];
  
V  = [-1.0, -1.0/sqrt(3.0), -1.0/sqrt(6.0); ...
       0.0,  2.0/sqrt(3.0), -1.0/sqrt(6.0); ...
       1.0, -1.0/sqrt(3.0), -1.0/sqrt(6.0); ...
       0.0,  0.0          ,  3.0/sqrt(6.0)];
   
vol = sqrt(sqrt(8.0)/3.0);

nK = G.cells.num;
nq = size(Xq,1);
      
I = zeros(nK,size(f([0,0,0]),2));

for i = 1:nK
    
    nodeNum = G.cells.nodePos(i):G.cells.nodePos(i+1)-1;
    nodes = G.cells.nodes(nodeNum);
    X = G.nodes.coords(nodes,:);
    tri = delaunay(X);
    nTri = size(tri,1);
    
    X = X(reshape(tri',[],1),:);
    ii = mcolon(1:4:4*nTri, 3:4:4*nTri);
    Xdiff = mat2cell(X(ii,:) - X(ii+1,:), 3*ones(nTri,1), 3);
    Xb = mat2cell(X(1:4:end,:),ones(nTri,1), 3);
    Vdiff = V(1:end-1,:) - V(2:end,:);
    
    A = cellfun(@(X) Vdiff\X, Xdiff, 'UniformOutput', false);
    b = cellfun(@(X,Y) X - V(1,:)*Y, Xb, A, 'UniformOutput', false);
    
    D = cellfun(@(X) abs(det(X)), A);
    
    Xhat = cell2mat(cellfun(@(X,y) bsxfun(@plus, Xq*X,y), A, b, ...
                                   'UniformOutput', false));
    
    I(i,:) = vol*repmat(w,1,nTri).*rldecode(D,4*ones(nTri,1),1)'*f(Xhat);

end

end