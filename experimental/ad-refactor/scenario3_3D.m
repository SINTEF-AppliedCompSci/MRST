function scenario3_3D(savename)

    gravity on;
    moduleCheck('ad-fi');
    moduleCheck('ad-refactor');

    %% Grid and rock parameters
    xnum      = 100;
    znum      = 60;%30;
    depth     = 700;
    dome_h    = 160;
    exponent  = 1.2;
    thickness = 100;
    G         = cartGrid([xnum, 1, znum], [10000, 3000, thickness]);
    G         = adjustVerticalCoords(G, depth, depth + thickness, znum, exponent, dome_h);
    G         = computeGeometry(G);
    cnum      = G.cells.num;
    rock.perm = repmat(1100*milli*darcy, cnum, 1);
    rock.poro = repmat(0.2             , cnum, 1);
    
    %% temperature related parameters
    tempS = 273.15 + 6;
    tgrad = 40;

    %% fluid
    EOS = CO2props('rho_big_trunc', []);

    fluid.pcGW    = @(sG) 0;                         % no capillary pressure - sharp interface
    fluid.relPerm = @(sG) deal(1-sG, sG);            % linear relperm
    fluid.rhoGS   = 1.977 * kilogram / meter^3;
    fluid.rhoWS   = 1100  * kilogram / meter^3;
    fluid.bW      = @(p,t) ones(numel(double(p)),1); % constant density
    fluid.bG      = @(p,t) EOS.rho(p, t) ./ fluid.rhoGS;
    fluid.muW     = @(p,t) 6.5e-4;                   % constant viscosity
    fluid.muG     = @(p,t) 5.36108e-5;               % constant viscosity

    %% Wells and schedule
    tnum     = 200; %60;
    inum     = 49;
    tot_time = 200 * year;
    rate     = 2e6 * kilo * kilogram / year / fluid.rhoGS;

    schedule = struct('W', verticalWell([], G, rock, 50, 1, [], ...
                                        'type'   , 'rate'  ,    ...
                                        'radius' , 0.3     ,    ...
                                        'comp_i' , [0 0 1] ,    ...
                                        'val'    , rate    ,    ...
                                        'name'   , 'I')    ,    ...
                      'step', struct('val', diff(linspace(0, tot_time, tnum+1)), ...
                                     'control', [ones(inum,1); zeros(tnum-inum,1)]));

    %% Initialize state
    hpress = @(z) (1 * atm + z .* fluid.rhoWS .* norm(gravity));

    state = struct('pressure', hpress(G.cells.centroids(:,3)), ...
                   's'       , repmat([1 0], G.cells.num, 1));

    %% Define boundary conditions
    bIx = @(grid, side) boundaryFaceIndices(grid, side, [], [], []);
    bc  = pside([], G, 'XMin', hpress(G.faces.centroids(bIx(G, 'XMin'), 3)));
    bc  = pside(bc, G, 'XMax', hpress(G.faces.centroids(bIx(G, 'XMax'), 3)));


    %% Define model
    model = twoPhaseGasWaterModel(G, rock, fluid, tempS, tgrad); 

    %% Run schedule
    [wellSols, states] = runScheduleRefactor(state, model, schedule, 'bc', bc);


    %% Compute caprock values and save results 

    for i = 1:numel(states)
        figure(1); clf;
        plotCellData(G, states{i}.s(:,2));
        view(0,0);

        figure(2); clf;
        plotCellData(G, states{i}.pressure);
        view(0,0);
        colorbar;

        pause(0.5);
    end
    
    save(savename, 'states', 'G');
end

function G = adjustVerticalCoords(G, zmin, zmax, zcellnum, exponent, dome_height)

    % adjust cell heights to make thinner cells on top
    cell_spacing = linspace(0, 1, zcellnum+1);
    cell_spacing = cell_spacing .^ exponent;
    cell_spacing = cell_spacing * (zmax-zmin);
    cell_spacing = reshape(repmat(cell_spacing, prod(G.cartDims(1:2)+1), 1), [],1);
    
    % Adjust geometry of grid to make a dome (in 2D)
    xnn = G.cartDims(1) + 1; % number of nodes along x direction
    dh   = dome_height / 2;  % half-height
    geom = dh * cos([1:xnn]'./xnn * 2 * pi);
    geom = (geom + fliplr(geom))/2;
    geom = reshape(repmat(geom, 1, prod(G.cartDims(2:3)+1)), [], 1);
    
    % set depth of dome ape
    z = geom + cell_spacing;
    z = z + (zmin-min(z));
    
    G.nodes.coords(:,3) = z;
end

        