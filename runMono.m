% No-flow boundary conditions, one source and one sink
close all;clear;clc

Lx=1;Ly=1;nx=11;ny=11;
G=cartGrid([nx ny],[Lx,Ly]);G=computeGeometry(G);
srcCells = findEnclosingCell(G,[3.5/11*Lx 0.5*Ly; 7.5/11*Lx 0.5*Ly]);
G=removeCells(G,srcCells);
figure, plotGrid(G);axis image
text(3.5/11*Lx-1/25,0.5*Ly,'\itp\rm=0');
text(7.5/11*Lx-1/25,0.5*Ly,'\itp\rm=1');

rock=makeUniformRock(G,[1000 1],67.5);
bc.face=boundaryFaces(G);bc.type=cell(numel(bc.face),1);
bc.value=cell(numel(bc.face),1);x=G.faces.centroids(bc.face,:);
d=bsxfun(@minus,x,[3.5/11*Lx 0.5*Ly]);d=sqrt(dot(d,d,2));ind1=d<Lx/nx;
d=bsxfun(@minus,x,[7.5/11*Lx 0.5*Ly]);d=sqrt(dot(d,d,2));ind2=d<Lx/nx;
bc.type(ind1)={'pressure'};bc.value(ind1)={@(x)0};
bc.type(ind2)={'pressure'};bc.value(ind2)={@(x)1};
bc.type(~ind1&~ind2)={'flux'};bc.value(~ind1&~ind2)={@(x)0};
src=addSource([],1,0);

s0=OnePhaseIncompTPFA(G,TransTPFA(G,rock,bc),src);
plotSurface(G.cells.centroids,s0.pressure,[nx ny],'griddata','TPFA');
s1=OnePhaseIncompMPFA(G,TransMPFAO(G,rock,bc),src);

interpFace=findHAP(G,rock,bc);
interpFace=correctHAP(G,interpFace);
OSflux=findOSflux(G,rock,bc,interpFace);
u0=ones(G.cells.num,1);
s2=PicardNTPFA(G,bc,src,OSflux,u0,1e-7,500);
s3=PicardNMPFA(G,bc,src,OSflux,u0,1e-7,500);

plotSurface(G.cells.centroids,s1.pressure,[nx ny],'griddata','MPFA-O');
plotSurface(G.cells.centroids,s2.pressure,[nx ny],'griddata','NTPFA');
plotSurface(G.cells.centroids,s3.pressure,[nx ny],'griddata','NMPFA');
disp(min([s1.pressure s2.pressure s3.pressure]))
disp(max([s1.pressure s2.pressure s3.pressure]))
disp([s2.iter s3.iter]);
