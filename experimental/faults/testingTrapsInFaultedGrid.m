%% Test trap analysis with faulted grid




%% Make grid and rock data
% We consider a sandbox of dimensions 10 km x 5 km x 50 m, that lies at a
% depth of 1000 m and has an inclination in the x-direction. For the top
% surface, we add a smooth sin/cos perturbation that will create domes. The
% porosity is set uniformly to 0.25
[Lx,Ly,H] = deal(10000, 5000, 50);
G  = cartGrid([20 20 1],[Lx Ly H]); % @@ originally 100x100x1
x  = G.nodes.coords(1:G.nodes.num/2,1)/Lx; 
y  = G.nodes.coords(1:G.nodes.num/2,2)/Ly;
z  = G.nodes.coords(1:G.nodes.num/2,3)/H;
zt = z + x - 0.2*sin(5*pi*x).*sin(5*pi*y.^1.5) - 0.075*sin(1.25*pi*y) + 0.15*sin(x+y);
zb = 1 + x;
G.nodes.coords(:,3) = [zt; zb]*H+1000;
G = computeGeometry(G);

%% Construct top surface grid
% We create a hybrid grid that represents the top surface. The grid is a 2D
% grid defined as a surface in 3D and contains a mapping from the 2D cells
% on the surfrace to the column of volumetric cells that lie below in the
% 3D grid. For visualization purposes, we create an extra grid that lies
% above the true top surface.
Gt = topSurfaceGrid(G);

figure;
Gt_zshifted = Gt; 
Gt_zshifted.nodes.z = Gt_zshifted.nodes.z - 15;
plot_opts = {'edgeColor', 'k', 'edgeAlpha', 0.1};
plotGrid(G, plot_opts{:});
plotCellData(Gt_zshifted, Gt_zshifted.cells.z, plot_opts{:});
view(30,25); axis tight

%% Geometric analysis of caprock (spill-point analysis)
% Compute traps and spill paths connecting them. Here, we use the
% cell-based algorithm. The cells that belong to the identified traps are
% colored white in the plot
res = trapAnalysis(Gt, true);
num_traps = max(res.traps);
plotGrid(Gt_zshifted, res.traps>0, 'FaceColor','white', 'EdgeAlpha',.1);


%% Show connection between traps and spill paths
% We make a 2D plot of the top surface in which traps are colored red,
% cells that lie along the connecting spill paths are colored green, and
% the remaining cells are colored blue. In addition, we display the number
% associated with each trap slightly above its local maximum.
clf
fpos = get(gcf,'Position');
set(gcf,'Position',[300 400 800 420],'PaperPositionMode','auto');
trap_field = zeros(size(res.traps));
trap_field(res.traps>0) = 2;
for r = [res.cell_lines{:}]'
    for c = 1:numel(r);
        trap_field(r{c}) = 1;
    end
end

subplot(2,3,[1 4]);
plotCellData(Gt, trap_field, 'EdgeColor','none');
view(90,90); axis tight off

for i=1:num_traps
   ind = res.top(i);
   text(Gt.cells.centroids(ind,1),Gt.cells.centroids(ind,2), res.trap_z(i)-50,...
      num2str(res.traps(ind)), 'Color',.99*[1 1 1], ...
      'FontSize', 14, 'HorizontalAlignment','center');
end
colormap(jet);



%% Introduce fault into grid
% specify the coordinates of the fault, convert it into face index
Xcoord = [1:Lx]';
Ycoord = 3000*ones(numel(Xcoord),1);


%[ ~, faultFaces ] = implementFaults(Gt, [], 'faultCoords',[Xcoords(:), Ycoords(:)]);




