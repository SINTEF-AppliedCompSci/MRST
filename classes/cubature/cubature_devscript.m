n = 10;
G = computeGeometry(cartGrid([n,n], [10,10]));

cubTri = TriangleCubature(G, 2);
cubLin = LineCubature(G, 2);