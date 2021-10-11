%%
dir='/home/hnil/Documents/data/benchmark_new'
ifile='ADGPRS_Final_Test/TNO_BENCHMARK_FINAL/TNO_BENCHMARK_FINAL.DATA'
deck_file=fullfile(dir,ifile)
mrstModule add deckformat
deck  = readEclipseDeck(deck_file);
%%
G=initEclipseGrid(deck)
clf,plotGrid(G)
%%
grd=deck.GRID;
grd_all=grd;
[i,j,k]=ind2sub(G.cartDims,G.cells.indexMap);
ijsub=unique([i,j],'rows');
grd_all.ACTNUM=ones(size(grd.ACTNUM));
ijind=sub2ind(G.cartDims(1:2),ijsub(:,1),ijsub(:,2));
ac=false(G.cartDims(1:2));
ac(ijind)=true;
ac=repmat(reshape(ac,[],1),G.cartDims(3),1);
grd_all.ACTNUM=ac;
G_a=processGRDECL(grd_all);
clf,plotGrid(G_a,'FaceColor','none')
plotGrid(G_a,find(grd.ACTNUM(G_a.cells.indexMap)==false))
%%

%%
mrstModule add project-mechanics-fractures
%%
grdc=grd_all;
grdc=padGrdecl(grdc,[0 0 1],[0 10000;0 5000;[1000 500]/50]);
grdc = verticalGrdecl(grdc);
%ac=false(grdc.cartDims(1:2));
%ac(ijind)=true;
%ac=repmat(reshape(ac,[],1),grdc.cartDims(3),1);
%grdc.ACTNUM=ac;
grdc = refineGrdeclLayers(grdc, [1 1], 10);
grdc = refineGrdeclLayers(grdc, [grdc.cartDims(3) grdc.cartDims(3)], 10);
ac=false(grdc.cartDims(1:2));
ac(ijind)=true;
ac=repmat(reshape(ac,[],1),grdc.cartDims(3),1);
grdc.ACTNUM=ac;
Gp=processGRDECL(grdc,'RepairZCORN',true);
clf,plotGrid(Gp)
%Gp = createAugmentedGrid(Gp);

%%
%grdc=cutGrdecl(grd_all,[45 118;50 138;1 16]);
grdc=cutGrdecl(grd_all,[45 118;50 138;1 16]);
%grdc=cutGrdecl(grd_all,[45 70;50 70;1 16]);
grdc=padGrdecl(grdc,[0 0 1],[0 10000;0 5000;[1000 500]/50]);
grdc = verticalGrdecl(grdc);
ref=4;
grdc = refineGrdeclLayers(grdc, [1 1], ref);
grdc = refineGrdeclLayers(grdc, [grdc.cartDims(3) grdc.cartDims(3)], ref);
G_cut=processGRDECL(grdc,'RepairZCORN',true,'Tolerance',1);
clf,plotGrid(G_cut),view(3)
plotGrid(G,'FaceColor','none')

mrstModule add vemmech
%%
profile off;profile on
G_cut=computeGeometry(G_cut);
Ga = createAugmentedGrid(G_cut);
profile off;profile viewer
%%
%disp_node = zeros(size(Ga.nodes.coords));
%disp_faces = zeros(size(Ga.faces.centroids(faces, :)));
% set all to free
%mask=false(size(disp_node));
%fix lower boundary
zmax = max(Ga.nodes.coords(:,3));
nodes = find(zmax==Ga.nodes.coords(:,3));
mask = true(numel(nodes),Ga.griddim);
disp_node =zeros(numel(nodes),Ga.griddim);
disp_bc = struct('nodes', nodes, 'uu', disp_node, 'faces', [], 'uu_face', [], 'mask', mask);
force_bc=[];
%% 
%
opt=struct('E',1e9,'nu',0.3);
Ev     = repmat(opt.E, Ga.cells.num, 1);
nuv    = repmat(opt.nu, Ga.cells.num, 1);
C      = Enu2C(Ev, nuv, Ga);

el_bc = struct('disp_bc', disp_bc, 'force_bc', force_bc);
lload=@(x) 0*x;
bbsize = 30000-(G.griddim-2)*20000;
lsolve = @mldivide;
profile off;profile on;
[uu,extra] = VEM_linElast(Ga, C, el_bc, lload, 'linsolve', lsolve, 'blocksize', bbsize, ...
                  'force_method', 'dual_grad_type','pressure',50*barsa*ones(Ga.cells.num,1),'experimental_scaling', true) ;
profile off;profile viewer

%%
figure(33)
clf,plotNodeDataDeformed(Ga,uu(:,3),uu*10),colorbar,view(33)

%%
A=extra.A;b=extra.rhs;
%profile off;profile on
disp('mldivide')
tic;
%linsolve=@mldivide;
%linsolve=@pcg(A,rhs,
x=A\b;
toc;
%{
disp('pcg ichol')
tic;
%L = ichol(A,struct('type','ict','droptol',1e-03,'michol','on'));
opts.type = 'nofill';
opts.michol = 'on';
L = ichol(A,opts);
[x1,flag1,rr1,iter1,rv1] = pcg(A,b,tol,maxit,L,L');
toc;
%}
%%
%%{
disp('bicg ilu')
tic;
setup.type = 'nofill';
setup.milu = 'row';
setup.droptol = 0.1;
[L,U] = ilu(A,setup);
toc;
[x1,flag1,rr1,iter1,rv1] = pcg(A,b,1e-5,40,L,U);
%}
%profile off;profile viewer

%%
% <html>
% <p><font size="-1">
% Copyright 2009-2021 SINTEF Digital, Mathematics & Cybernetics.
% </font></p>
% <p><font size="-1">
% This file is part of The MATLAB Reservoir Simulation Toolbox (MRST).
% </font></p>
% <p><font size="-1">
% MRST is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
% </font></p>
% <p><font size="-1">
% MRST is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
% </font></p>
% <p><font size="-1">
% You should have received a copy of the GNU General Public License
% along with MRST.  If not, see
% <a href="http://www.gnu.org/licenses/">http://www.gnu.org/licenses</a>.
% </font></p>
% </html>
