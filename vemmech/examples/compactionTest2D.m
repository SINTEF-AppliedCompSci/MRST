
%% Define parameters

opt = struct('grid_type', 'cartgrid', ... % see squareGrid
            'disturb'   , 0.0, ...      % parameter for disturbing grid
            'E'         , 0.3*1e9, ...  % young's modulus
            'nu'        , 0.4, ...
            'islinear'  , false);       %poisson's ratio

opt.cartDims       = [10 10];
opt.L              = [15*10/10 15];
opt.hanging        = false;
opt.free_top       = true;
opt.triangulate    = true; % triangulate some faces, see below
opt.gravity_load   = true;
opt.top_load       = false;
opt.gtol           = 0.1e-1;
% opt.ref            = 10;
opt.twist          = true;
opt.disturb        = 0.05;
% Two different methods are implemented to compute the loading term, see paper
opt.force_method   = 'cell_force'; % 'dual_grad_type'

%% Construct grid

G = squareGrid(opt.cartDims, opt.L, 'grid_type', opt.grid_type);

if(opt.triangulate)
    face = unique(G.cells.faces(any(bsxfun(@eq, G.cells.faces(:, 2), [3, 4]), 2), 1));
    G = triangulateFaces(G, face);
    G = sortEdges(G);
end

if(opt.twist)
    G = twister(G, opt.disturb);
end

G = createAugmentedGrid(G);
G = computeGeometry(G);

figure()
clf, plotGrid(G, 'Marker', '*')

%% Setup loads

[el_bc, load] = makeCompactionTest(G, opt);


%% Define rock parameters

Ev  = repmat(opt.E, G.cells.num, 1);
nuv = repmat(opt.nu, G.cells.num, 1);
C   = Enu2C(Ev, nuv, G);

lsolve = @mldivide;

%% Assemble and solve the system

bbsize = 30000-(G.griddim-2)*20000; % block size for assembly
uu = VEM_linElast(G, C, el_bc, load, ...
                  'linsolve'    , lsolve, ...
                  'blocksize'   , bbsize, ...
                  'force_method', opt.force_method);

%% Assemble divergence operator, see paper

div = VEM_div(G);

%% Plotting

figure(), 
clf
plotNodeDataDeformed(G, uu(:, G.griddim), uu);
% cdiv = div*reshape(uu', [], 1);
% plotCellDataDeformed(G, cdiv, fac*uu, 'EdgeColor', 'none');colorbar
% plotFaces(G, bc{1}.face, 'FaceColor', 'r')
% plotFaces(G, bc{2}.face, 'FaceColor', 'b')
% set(cb, 'YTick', [-2e-2:1e-2:0])
% axis off
%t itle(mycase)

%% Compute  the analytical solution

if(~isempty(el_bc.force_bc))
    ff = abs(el_bc.force_bc.force(1, G.griddim));
else
    ff = 0;
end
start = max(G.faces.centroids(:, G.griddim));
top   = min(G.faces.centroids(:, G.griddim));

fac = 100*300/2;

ana = @(z) ff*(z-start)./(C(1, 1))-double(opt.gravity_load)*fac*((top-start).^2 - (z-top).^2)/C(1, 1);
divana = @(z) (ff./C(1, 1))-double(opt.gravity_load)*fac*(-2*(z-top))/C(1, 1);

%% Comparison plots

z = G.nodes.coords(:, G.griddim);
z(abs(ana(z))<max(abs(ana(z)))*1e-2) = nan;
zl = unique(z);
figure(), 
subplot(4, 1, 1)
plot(z, uu(:, G.griddim), '*', zl, ana(zl))
subplot(4, 1, 2)
err = (uu(:, G.griddim)-ana(z))./max(abs(ana(z)));
plot(z, err, '*')
subplot(4, 1, 3)
div = VEM_div(G);
plot(G.cells.centroids(:, G.griddim), div*reshape(uu', [], 1)./G.cells.volumes, '*');
subplot(4, 1, 4)
div = VEM_div(G);
zc = G.cells.centroids(:, G.griddim);
diverr = (div*reshape(uu', [], 1)./G.cells.volumes-divana(zc))./max(abs(divana(zc)));
plot(zc, diverr, '*');

