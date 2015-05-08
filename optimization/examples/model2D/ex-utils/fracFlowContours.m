function [] = fracFlowContours(G, W, states, fluid, frac, varargin)
% utility-function for analyseModel2D.m
% plot contour of each state where fracFlow = frac
% only for 2D rectangular grids and oil/water systems!
assert(or(G.griddim == 2, G.cartDims(3)==1));
assert(prod(G.cartDims)==G.cells.num);

mn = min(G.nodes.coords(:,1:2));
mx = max(G.nodes.coords(:,1:2));
dd = (mx-mn)./G.cartDims(1:2);
x = (mn(1):dd(1):mx(1))';
x = .5*(x(1:end-1)+x(2:end));
y = (mn(2):dd(2):mx(2))';
y = .5*(y(1:end-1)+y(2:end));

ns = numel(states);

for k = 1:ns
    set(gca,'FontSize', 14)
    hold on
    plotCellData(G, states{k}.s(:,2))
    s = states{k}.s(:,1);
    [krw, kro] = fluid.relPerm(s, 0*s);
    %[krw, kro] = fluid.relPerm(states{k}.s(:,1));
    muw = fluid.muW(states{k}.pressure);
    muo = fluid.muO(states{k}.pressure);
    bw = fluid.bW(200*barsa);   % FVF at reference pressure
    bo = fluid.bO(200*barsa); 
    fw = (bw*krw./muw)./(bw*krw./muw+bo*kro./muo);
    m = reshape(fw/frac, G.cartDims(1:2))';
    [c,h] = contour(x,y,m, [1 1], varargin{:});
    caxis([0 1])
    axis off, axis equal
    view([-1 -2 1.5])
    camproj perspective
    axis tight
    plotWell(G, W, 'Fontsize', 14, 'Color', 'k')
end
end


    



