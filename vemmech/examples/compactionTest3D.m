
%% Define parameters

griddim = 3;
opt = struct('L'         , [1000 1000 100], ...
             'cartDims'  , ones(1, griddim)*ceil((1e3).^(1/griddim)), ...
             'grid_type' , 'triangle', ...
             'disturb'   , 0.0, ...     % parameter for  grid distortion
             'E'         , 0.3*1e9, ... % Young's modulus
             'nu'        , 0.4, ...
             'islinear'  , false);      % Poison's ratio
opt.L            = [15 15 3];
opt.islinear     = false;
opt.force_method = 'dual_grad_type';
opt.hanging      = false;
opt.free_top     = true;
opt.use_pressure = false;
opt.triangulate  = true;% triangulate some faces
opt.vertical     = false;%only valid for norne
opt.gravity_load = true;
opt.top_load     = true;
opt.gtol         = 0.1e-1;
opt.ref          = 10;
grid_case        = 'box';
opt.cartDims     = [[1 1]*3 10];
opt.flipgrid     = false;

grid_case = 'box'; % Other choices: 'grdecl' 'sbed' 'norne' 'model3'

%% Construct grid

G = complex3DGrid(opt, grid_case);
if(opt.flipgrid)
    G = flipGrid(G);
end
G = mrstGridWithFullMappings(G);
G = computeGeometry(G);

figure()
clf;
plotGrid(G);


%% Setup loads

[el_bc, load] = makeCompactionTest(G, opt);


%% Define rock parameters
Ev = repmat(opt.E, G.cells.num, 1);
nuv = repmat(opt.nu, G.cells.num, 1);
C = Enu2C(Ev, nuv, G);
lsolve = @mldivide;

G = computeGeometry(G);

% Change load for pressure since it is a constant vector
if(opt.use_pressure)
    minp           = -max(sum(load(G.faces.centroids).*G.faces.centroids, 2));
    pressure       = @(X)(-sum(load(X).*X, 2)-minp);
    load           = @(X)(load(X)*0);
    bccoord        = G.faces.centroids(el_bc.force_bc.faces, :);
    press_boundary = sum(pressure(bccoord));
    bsign          = (2* (G.faces.neighbors(el_bc.force_bc.faces, 2) == 0)-1);
    normals        = -bsxfun(@times, G.faces.normals(el_bc.force_bc.faces, :), bsign);
    el_bc.force_bc.force = el_bc.force_bc.force+bsxfun(@times, normals, press_boundary);
else
    pressure = @(X)(-sum(load(X).*X, 2)*0.0);
end

%% Assemble and solve the system

bbsize = 30000-(G.griddim-2)*20000;
uu = VEM_linElast(G, C, el_bc, load, 'linsolve', lsolve, 'blocksize', bbsize, ...
                  'force_method', opt.force_method, 'pressure', ...
                  pressure(G.cells.centroids));

%% Assemble divergence operator, see paper

div = VEM_div(G);

%% Plotting

figure();
clf;
plotNodeDataDeformed(G, uu(:, 3), uu);cb = colorbar;view(3)
cdiv = div*reshape(uu', [], 1);
axis off

%% Compute  the analytical solution

ff = abs(el_bc.force_bc.force(1, 3));
start = max(G.faces.centroids(:, 3));
top = min(G.faces.centroids(:, 3));
[lambda, mu] = ENu2LMu_3D(opt.E, opt.nu);
ana = @(z) ff*(z-start)./(C(1, 1))+double(opt.gravity_load)*10*300*((z).^2-(start).^2)/C(1, 1);
ana = @(z) ff*(z-start)./(C(1, 1))-double(opt.gravity_load)*50*300*((top-start).^2 - (z-top).^2)/C(1, 1);
divana = @(z) (ff./C(1, 1))-double(opt.gravity_load)*50*300*(-2*(z-top))/C(1, 1);

%% Comparison plots

z = G.nodes.coords(:, 3);
z(abs(ana(z))<max(abs(ana(z)))*1e-2) = nan;
zl = unique(z);
figure(),
subplot(4, 1, 1)
plot(z, uu(:, 3), '*', zl, ana(zl))
subplot(4, 1, 2)
err = (uu(:, 3)-ana(z))./abs(ana(z));
plot(z, err, '*')
subplot(4, 1, 3)
div = VEM_div(G);
plot(G.cells.centroids(:, 3), div*reshape(uu', [], 1)./G.cells.volumes, '*');
subplot(4, 1, 4)
div = VEM_div(G);
zc = G.cells.centroids(:, 3);
diverr = (div*reshape(uu', [], 1)./G.cells.volumes-divana(zc))./abs(divana(zc));
plot(zc, diverr, '*');
