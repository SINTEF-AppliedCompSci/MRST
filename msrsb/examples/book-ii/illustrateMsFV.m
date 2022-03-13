%% Illustration of the MsFV method
% This example illustrates the algebraic construction of the multiscale
% finite-volume (MsFV) method, starting with the wirebasket ordering and
% the structure of the fine-scale matrix before and after permutation. We
% then show how to construct the prolongation and restriction operators and
% show their sparsity pattern along with that of the reduced coarse-scale
% system.

mrstModule add coarsegrid incomp spe10

%% Set up model
% The mesh is designed to have k*k coarse blocks, that each consists of
% (2m+1)*(2m+1) cells with petrophysical properties sampled from the 15th
% layer of Model 2 of the 10th SPE Comparative Solution Project. By setting
% the flag 'MatrixOutput' to true, we get the linear system as output in
% the fields state.A and state.rhs.
k = 3;
m = 2;
n = k*(2*m+1);
G      = computeGeometry(cartGrid([n n]));
rock   = getSPE10rock(1:n, 1:n, 15);
rock.perm = rock.perm(:, 1);
hT     = computeTrans(G, rock);
bc     = pside([], G, 'West', 100*barsa);
bc     = pside(bc, G, 'East',  50*barsa);
fluid  = initSingleFluid('rho', 1, 'mu', 1);
state0 = initResSol(G, 0);
state  = incompTPFA(state0, G, hT, fluid, 'MatrixOutput', true, 'bc', bc);

%% Show partition and wirebasket ordering
% For this simple and symmetric geometry, the partition can be computed
% explicitly by marking the appropriate edge rows in an indicator matrix D
% and then adding D and D'. The resulting matrix will have value 2 on
% nodes, 1 on edges and zero in inner cells.
p = partitionUI(G, [k,k]);
d = zeros(n);
d(:,m+1:2*m+1:n) = 1;
d = d + d';
figure(1), clf, plotCellData(G,d(:),'EdgeAlpha',.1); axis tight
colormap(.6*jet+.4*ones(size(jet)));
outlineCoarseGrid(G, p,'Color','k','LineWidth',1); axis off
title('Wirebasket ordering')

%% Alternative visualization of wirebasket ordering
% Show the porosity and outline edge and node cells.
clf,
h1=plotCellData(G,rock.poro,'EdgeAlpha',.1); axis tight
outlineCoarseGrid(G, p,'Color','k','LineWidth',2);
h2=plotGrid(G,d(:)==1,'EdgeColor','w','FaceColor','w','LineWidth',2,'FaceAlpha',0);
h3=plotGrid(G,d(:)==2,'EdgeColor','r','FaceColor','w','LineWidth',2,'FaceAlpha',0);
colormap(parula); axis off
h=legend([h3,h2,h1],' node',' edge',' inner');
set(h,'Color',[.8 .8 .9],'FontSize',14);
h.Position = [.65 .67 .25 .25];
title('Wirebasket ordering')

%% Show matrix and wirebasket reordering
% We show the sparsity pattern of the fine-scale matrix before and after
% wirebasket reordering
A = state.A;
i = find(d(:)==0)'; in=numel(i);
e = find(d(:)==1)'; en=numel(e);
n = find(d(:)==2)'; nn=numel(n);
figure(2), clf
subplot(1,2,1), spy(A); title('Original matrix A');
subplot(1,2,2), cla, hold on
patch([0 in in 0], [0 0 in in],[0 0 0 0],'FaceColor',[.8 1 .8],'Edgecolor','none');
patch([0 en en 0]+in, [0 0 en en]+in,[0 0 0 0],'FaceColor',[.8 .8 1],'Edgecolor','none');
patch([0 in in 0], [0 0 en en]+in,[0 0 0 0],'FaceColor',[.8 1 1],'Edgecolor','none');
patch([0 en en 0]+in, [0 0 in in],[0 0 0 0],'FaceColor',[.8 1 1],'Edgecolor','none');
patch([0 nn nn 0]+in+en, [0 0 nn nn]+in+en,[0 0 0 0],'FaceColor',[1 .6 .6],'Edgecolor','none');
patch([0 in+en in+en 0], [0 0 nn nn]+in+en,[0 0 0 0],'FaceColor',[1 1 .8],'Edgecolor','none');
patch([0 0 nn nn]+in+en, [0 in+en in+en 0], [0 0 0 0],'FaceColor',[1 1 .8],'Edgecolor','none');
patch([in in+en in+en in], [0 0 nn nn]+in+en, [0 0 0 0],'FaceColor',[1 .8 1],'Edgecolor','none');
patch([0 0 nn nn]+in+en, [in in+en in+en in], [0 0 0 0],'FaceColor',[1 .8 1],'Edgecolor','none');
spy(A([i e n],[i e n]));
title('After wirebasket ordering');

%% Construct and visualize P
% Using the explicit formulas given in the book chapter, we construct the
% prolongation operator and show its sparsity structure in wirebasket and
% natural cell ordering. We also show one of the basis functions on top of
% the 2D mesh.
Aii = A(i,i);
Aie = A(i,e);
Aen = A(e,n);
Mee = A(e,e)+diag(Aie');
Pn = diag(ones(nn,1));
Pe = -inv(Mee)*Aen;
Pi  = Aii\(Aie*(-Pe));
Pwb = full([Pi; Pe; Pn]');

% Show in wirebasket ordering
figure, subplot(1,2,1)
g = cartGrid(size(Pwb)); val=Pwb(:); val(val==0)=NaN;
plotCellData(g,val,'EdgeColor','none'); 
axis equal tight, box on, caxis([0 1]);
set(gca,'YDir','reverse','XTick',[],'YTick',[],'DataAspectRatio',[1 3 1]);
P = Pwb; P(:,[i e n])=Pwb;
title('P, wirebasket ordering');

% Show in natural ordering
subplot(1,2,2); val=P(:); val(val==0)=NaN;
plotCellData(g, val, 'EdgeColor','none'); 
axis equal tight, box on, caxis([0 1]);
set(gca,'YDir','reverse','XTick',[],'YTick',[],'DataAspectRatio',[1 3 1]);
title('P, natural ordering');

% Show on top of grid
figure, plotCellData(G,rock.poro,'EdgeAlpha',.1);
outlineCoarseGrid(G, p,'Color','k','LineWidth',1);
plotCellData(G,max(rock.poro)*P(5,:)',P(5,:)'>0,'EdgeAlpha',.1);
title('Basis function imposed on porosity');

%% Construct and visualize R
% The restriction operator can either be set as P' or as a finite-volume
% summation operator. Show the sparsity structure of both.
figure
subplot(2,1,1);
R = P';
g = cartGrid(size(R)); val=R(:); val(val==0)=NaN;
plotCellData(g, val, 'EdgeColor','none'); 
axis equal tight, box on, caxis([0 1]);
set(gca,'YDir','reverse','XTick',[],'YTick',[],'DataAspectRatio',[3 1 1]);
title('Galerkin restriction operator');

subplot(2,1,2)
pval=unique(p);
p=p';
for i=1:numel(pval)
    R(:,i) = (p==pval(i))+0;
end
val=R(:); val(val==0)=NaN;
plotCellData(g, val, 'EdgeColor','none'); 
axis equal tight, box on, caxis([0 1]);
set(gca,'YDir','reverse','XTick',[],'YTick',[],'DataAspectRatio',[3  1 1]);
title('Finite-volume restriction operator');

%% Compare linear systems
% Show sparsity plot of the fine-scale and the coarse-scale systems
A = full(A); Ac = R'*A*P';
figure, 
subplot(1,2,1),val=A(:); val(val==0)=NaN;
plotCellData(cartGrid(size(A)), val, 'EdgeColor','none');
axis equal tight, box on, set(gca,'YDir','reverse','XTick',[],'YTick',[]);
title('Fine-scale discretization')

subplot(1,2,2)
val=Ac(:); val(val==0)=NaN;
plotCellData(cartGrid(size(Ac)), val, 'EdgeColor','none');
axis equal tight, box on, set(gca,'YDir','reverse','XTick',[],'YTick',[]);
title('Coarse-scale discretization');
%%
% <html>
% <p><font size="-1">
% Copyright 2009-2022 SINTEF Digital, Mathematics & Cybernetics.
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
