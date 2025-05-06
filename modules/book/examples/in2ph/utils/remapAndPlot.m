function [x] = remapAndPlot(Gi, xi, mx, cval, varargin)
% Remap between a parallel and a diagonal grid and plot contours
%
if nargin>=5
  fig = gcf;
  figure
  subplot(1,2,1),
  plotCellData(Gi,xi), axis equal tight;
  hold on; 
  plot([0 0 mx([1 1]) 0],[0 mx([2 2]) 0 0],'k','LineWidth',2); hold off
end

theta = pi/4;
R = [cos(theta) -sin(theta); sin(theta) cos(theta)];
G = cartGrid(Gi.cartDims.*[1 2 1], [mx 1].*[1 2 1]);
G.nodes.coords(:,2) = G.nodes.coords(:,2) - mx(2);
G.nodes.coords(:,1:2) = sqrt(2)*(R*(G.nodes.coords(:,1:2)'))';
G = computeGeometry(G);
x = reshape(xi, Gi.cartDims);
x = x(:,[end:-1:1 1:end],:);
x = x(:);

if nargin>=5
    subplot(1,2,2);
    plotCellData(G,x); axis equal tight
    hold on; 
    plot([0 0 mx([1 1]) 0],[0 mx([2 2]) 0 0],'r','LineWidth',2); hold off
    figure(fig);
end

contourf(reshape(G.cells.centroids(:,1), G.cartDims),...
    reshape(G.cells.centroids(:,2), G.cartDims), ...
    reshape(x(:),G.cartDims), [0 cval 1],'EdgeColor','none');
% plotCellData(G,x,'EdgeColor','none'); caxis([-.05 1.05]);
end

