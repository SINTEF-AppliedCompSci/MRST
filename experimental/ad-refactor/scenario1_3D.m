function scenario1_3D(savename)

    gravity on;
    moduleCheck('ad-fi');
    moduleCheck('ad-refactor');

    %% Grid and rock parameters
    znum      = 30;%60;%30;
    depth     = 750;
    thickness = 150;
    exponent  = 1.2;
    G         = cartGrid([25, 1, znum], [40000, 3000, thickness]);
    G         = adjustVerticalCoords(G, depth, depth + thickness, znum, exponent);
    G         = computeGeometry(G);
    cnum      = G.cells.num;
    rock.perm = repmat(400*milli*darcy, cnum, 1);
    rock.poro = repmat(0.1            , cnum, 1);
    
    %% temperature related parameters
    tempS = 273.15 + 6;
    tgrad = 40;

    %% fluid
    EOS = CO2props('rho_big_trunc', []);

    fluid.pcGW    = @(sG) 0;                         % no capillary pressure - sharp interface
    fluid.relPerm = @(sG) deal(1-sG, sG);            % linear relperm
    fluid.rhoGS   = 1.977 * kilogram / meter^3;
    fluid.rhoWS   = 1050  * kilogram / meter^3;
    fluid.bW      = @(p,t) ones(numel(double(p)),1); % constant density
    fluid.bG      = @(p,t) EOS.rho(p, t) ./ fluid.rhoGS;
    fluid.muW     = @(p,t) 5.4e-5;                   % constant viscosity
    fluid.muG     = @(p,t) 5.36108e-5;               % constant viscosity

    %% Wells and schedule
    tnum     = 60; %60;
    inum     = 19;
    tot_time = 60 * year;
    rate     = 1e7 * kilo * kilogram / year / fluid.rhoGS;
    wcell    = ceil(G.cartDims(1:2)/2);

    schedule = struct('W', verticalWell([], G, rock, wcell(1), wcell(2), [], ...
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

    % for i = 1:numel(states)
    %     figure(1); clf;
    %     plotCellData(G, states{i}.s(:,2));
    %     view(0,0);

    %     figure(2); clf;
    %     plotCellData(G, states{i}.pressure);
    %     view(0,0);
    %     colorbar;

    %     pause(0.5);
    % end

    result = makeResultStructure(states, G, rock, EOS, fluid.rhoWS, @(d) tempS + d * tgrad/1000);
    save(savename, 'result', 'G');
end
% ----------------------------------------------------------------------------
function G = adjustVerticalCoords(G, zmin, zmax, zcellnum, exponent)
    
    tmp = linspace(0, 1, zcellnum+1);
    tmp = tmp .^ exponent;
    tmp = tmp * (zmax-zmin) + zmin;
    tmp = reshape(repmat(tmp, prod(G.cartDims(1:2)+1), 1), [],1);
    G.nodes.coords(:,3) = tmp;
end

        