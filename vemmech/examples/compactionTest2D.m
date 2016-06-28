
%% Define parameters

griddim = 3;
opt = struct('L'          , [1000 1000 100], ...
            'cartDims'  , ones(1, griddim)*ceil((1e3).^(1/griddim)), ...
            'grid_type' , 'triangle', ...
            'disturb'   , 0.0, ...      % parameter for disturbing grid
            'E'         , 0.3*1e9, ...  % young's modulus
            'nu'        , 0.4, ...
            'islinear'  , false);       %poisson's ratio

opt.L              = [15 15 3];
opt.hanging        = false;
opt.free_top       = true;
opt.triangulate    = true; % triangulate some faces
opt.vertical       = false;   % only valid for norne
opt.gravity_load   = true;
opt.top_load       = false;
opt.gtol           = 0.1e-1;
opt.ref            = 10;
opt.pressure       = false;
opt.mycase         = 'cartgrid';
opt.twist          = true;opt.disturb = 0.05;
% Two different methods are implemented to compute the loading term, see paper
opt.force_method   = 'cell_force'; % 'dual_grad_type'
opt.grid_type      = 'cartgrid';
opt.cartDims       = [10 10];
opt.L              = [15*10/10 15];
opt.flipgrid       = false;

%% Construct grid

G = squareGrid(opt.cartDims, opt.L, 'grid_type', opt.grid_type);

if(opt.triangulate)
    face = unique(G.cells.faces(any(bsxfun(@eq, G.cells.faces(:, 2), [3, 4]), 2), 1));
    face = unique(G.cells.faces(any(bsxfun(@eq, G.cells.faces(:, 2), [1, 4]), 2), 1));
    G = triangulateFaces(G, face);
    G = sortEdges(G);
end

if(opt.twist)
    G = twister(G, opt.disturb);
end

if(opt.flipgrid)
    G = flipGrid(G);
end

G = mrstGridWithFullMappings(G);
G = computeGeometryCalc(G);

figure()
clf, plotGrid(G, 'Marker', '*')

%% Setup loads

[el_bc, load] = makeCompactionTest(G, opt);


%% Define rock parameters

Ev  = repmat(opt.E, G.cells.num, 1);
nuv = repmat(opt.nu, G.cells.num, 1);
C   = Enu2C(Ev, nuv, G);

lsolve = @mldivide;
if(strcmp(mycase, 'box'))
    lsolve = @agmg;
end

G = computeGeometry(G);

% change load for pressure since it is a constant vector
minp = -max(sum(load(G.faces.centroids).*G.faces.centroids, 2));

if(opt.pressure)
    pressure = @(X)(-sum(load(X).*X, 2) - minp);
    load = @(X) load(X)*0;
    %{
    bccoord = G.faces.centroids(el_bc.force_bc.faces, :);
    press_boundary = sum(pressure(bccoord));
    bsign = (2* (G.faces.neighbors(el_bc.force_bc.faces, 2) == 0)-1);
    normals = -bsxfun(@times, G.faces.normals(el_bc.force_bc.faces, :), bsign);
    el_bc.force_bc.force = el_bc.force_bc.force+bsxfun(@times, normals, press_boundary);
    %}
else
    pressure = @(X)(-sum(load(X).*X, 2)*0.0);  
end

%% Assemble and solve the system

bbsize = 30000-(G.griddim-2)*20000; % block size for assembly
uu = VEM_linElast(G, C, el_bc, load, 'linsolve', lsolve, 'blocksize', bbsize, 'force_method', opt.force_method, 'pressure', pressure(G.cells.centroids));

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
top = min(G.faces.centroids(:, G.griddim));

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

