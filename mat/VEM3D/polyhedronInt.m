function I = polyhedronInt(G,cells,f, k)

%   http://people.sc.fsu.edu/~jburkardt/m_src/tetrahedron_arbq_rule/tetrahedron_arbq_rule.html

if k == 2

Xq = [-.2677028585910073,  .5939017006199488E-01, -.3969426941142150; ...
       .1510931841533142,  .4354087732476309,      .2151719254306919; ...
      -.1367699195237390, -.3280341115590410,      .2950846519937133; ...
       .8067449309051964, -.3400977314285288,     -.3430378951002550];

w  = [ 0.2918008865477151, 0.2706884392356724, ...
       0.3098349601753470, 0.9865925745591227E-01 ];
   
else
    xs = [ ...
    -.1685037180276000,0.2783799427534418E-01, ...
    -.3512216177343445,0.4308532463549043, ...
    -.4676763747967377,0.1492831253848431 ];
  ys = [ ...
    0.1910914916271708,-.2304932838839657E-01, ...
    0.1835144026339993,-.2474715823180446, ...
    -.4235250827264375,0.6397847685164516 ];
  zs = [ ...
    -.3896267314585163,0.5481350663241830, ...
    0.5147815330343534E-01,-.1683315641007033, ...
    -.1586973077889307,-.1080219253055393 ];
  w = [ ...
    0.1287213727402025,0.2179034339695993, ...
    0.1243503792062836,0.2446917182410072, ...
    0.1365439875826512,0.1187726516749031 ];
    Xq = [xs(:), ys(:), zs(:)];
end
    
  
V  = [-1.0, -1.0/sqrt(3.0), -1.0/sqrt(6.0); ...
       0.0,  2.0/sqrt(3.0), -1.0/sqrt(6.0); ...
       1.0, -1.0/sqrt(3.0), -1.0/sqrt(6.0); ...
       0.0,  0.0          ,  3.0/sqrt(6.0)];   
Vdiff = V(1:end-1,:) - V(2:end,:);
   
vol = sqrt(sqrt(8.0)/3.0);

nK = numel(cells);
nq = size(Xq,1);
      
I = zeros(nK,size(f([0,0,0]),2));

for i = 1:nK
    
    nodeNum = G.cells.nodePos(cells(i)):G.cells.nodePos(cells(i)+1)-1;
    nodes = G.cells.nodes(nodeNum);
    X = G.nodes.coords(nodes,:);
    tri = delaunay(X);
    nTri = size(tri,1);
    
    X = X(reshape(tri',[],1),:);
    ii = mcolon(1:4:4*nTri, 3:4:4*nTri);
    Xdiff = mat2cell(X(ii,:) - X(ii+1,:), 3*ones(nTri,1), 3);
    Xb = mat2cell(X(1:4:end,:),ones(nTri,1), 3);

    A = cellfun(@(X) Vdiff\X, Xdiff, 'UniformOutput', false);
    b = cellfun(@(X,Y) X - V(1,:)*Y, Xb, A, 'UniformOutput', false);
    
    D = cellfun(@(X) abs(det(X)), A);
    
    Xhat = cell2mat(cellfun(@(X,y) bsxfun(@plus, Xq*X,y), A, b, ...
                                   'UniformOutput', false));
    
    I(i,:) = vol*repmat(w,1,nTri).*rldecode(D,nq*ones(nTri,1),1)'*f(Xhat);

end

end