function plotObsAndSim(Gt, obs, sim, varargin)
% Plot simulated and observed plume thicknesses
% obs and sim are cell arrays containing plume thicknesses (heights) "h"

    assert(numel(sim) == 12) % 12 for 1999-2010
    assert(numel(obs) == 12)

    opt.edgealpha = 0.1;
    opt.edgecolor = 'none';
    opt = merge_options(opt, varargin{:});
    
    myplotGrid = @(G) plotGrid(G, 'facecolor','none','edgealpha',opt.edgealpha);
    myplotCellData = @(G,data,tol) plotCellData(G, data, data>tol, ...
        'edgealpha',opt.edgealpha, 'edgecolor',opt.edgecolor);

    %ty=[3,6,8,12];
    Tol = 0.01;

    figure,
    subplot(2,4,1)
    title({'Year: 2001';'';'simulated'})
    myplotGrid(Gt);
    myplotCellData(Gt, sim{3}.h, Tol); colormap(jet); axis equal tight off; colorbar
    subplot(2,4,5)
    title('observed')
    myplotGrid(Gt);
    myplotCellData(Gt, obs{3}.h, Tol); colormap(jet); axis equal tight off; colorbar

    subplot(2,4,2)
    title({'Year: 2004';'';'simulated'})
    myplotGrid(Gt);
    myplotCellData(Gt, sim{6}.h, Tol); colormap(jet); axis equal tight off; colorbar
    subplot(2,4,6)
    title('observed')
    myplotGrid(Gt);
    myplotCellData(Gt, obs{6}.h, Tol); colormap(jet); axis equal tight off; colorbar

    subplot(2,4,3)
    title({'Year: 2006';'';'simulated'})
    myplotGrid(Gt);
    myplotCellData(Gt, sim{8}.h, Tol); colormap(jet); axis equal tight off; colorbar
    subplot(2,4,7)
    title('observed')
    myplotGrid(Gt);
    myplotCellData(Gt, obs{8}.h, Tol); colormap(jet); axis equal tight off; colorbar

    subplot(2,4,4)
    title({'Year: 2010';'';'simulated'})
    myplotGrid(Gt);
    myplotCellData(Gt, sim{12}.h, Tol); colormap(jet); axis equal tight off; colorbar
    subplot(2,4,8)
    title('observed')
    myplotGrid(Gt);
    myplotCellData(Gt, obs{12}.h, Tol); colormap(jet); axis equal tight off; colorbar

end
