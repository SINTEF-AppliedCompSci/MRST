function scenario2_3D(savename)

    gravity on;
    moduleCheck('ad-fi');
    moduleCheck('ad-refactor');

    %% Grid and rock parameters
    znum      = 60; %30
    depth     = 1000;
    thickness = 150;
    slope     = 1/2 * pi / 180;
    exponent  = 1;%1.3;
    G         = cartGrid([200, 1, znum], [80000, 3000, thickness]);
    G         = adjustVerticalCoords(G, depth, depth + thickness, znum, exponent);
    G         = rotateGrid(G, slope, G.nodes.coords(1,:)', [0 1 0]');
    G         = computeGeometry(G);
    cnum      = G.cells.num;
    rock.perm = repmat(1400*milli*darcy, cnum, 1);
    rock.poro = repmat(0.1             , cnum, 1);
    
    %% temperature related parameters
    tempS = 273.15 + 8;
    tgrad = 45;
    tinfo     = {tempS, 0, tgrad}; 

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
    fluid.sw      = 0; % residual water (required by the below call to 'height2Sat')
    fluid.sr      = 0; % residual CO2   (ditto)

    %% Wells and schedule
    tnum     = 50;
    tot_time = 50 * year;
    schedule = struct('W', [], ...
                      'step', struct('val', diff(linspace(0, tot_time, tnum+1)), ...
                                     'control', [zeros(tnum,1)]));

    %% Initialize state
    mass_profile = zeros(200,1);
    mass_profile(20:40) = linspace(0, 100, 21) * 4e7;
    mass_profile(40:60) = linspace(100, 0, 21) * 4e7;
    
    [s, p] = computeSaturationAndPressure(G, mass_profile, EOS, tinfo, fluid, ...
                                          rock.poro, 1*atm);

    state = struct('pressure', p, 's', s);
    
    %% Define boundary conditions
    hpress = @(z) (1 * atm + z .* fluid.rhoWS .* norm(gravity));
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

    %     % figure(2); clf;
    %     % plotCellData(G, states{i}.pressure);
    %     % view(0,0);
    %     % colorbar;

    %     pause(0.5);
    % end
    result = makeResultStructure(states, G, rock, EOS, fluid.rhoWS, @(d) tempS + d * tgrad/1000);       save(savename, 'result', 'G');
    %save(savename, 'states', 'G'); 
end

function [s, p] = computeSaturationAndPressure(G, mass_profile, EOS, tinfo, ...
                                               fluid, poro, top_press)
% Compute 3D saturations from a 2D mass profile
    
    %% Computing saturation of each cell
    Gt       = topSurfaceGrid(G);
    topPress = fluid.rhoWS * norm(gravity) * Gt.cells.z + 1*atm;
    tmp.h    = convertMassesToHeight(mass_profile, EOS, topPress, tinfo, fluid.rhoWS, ...
                                     0, Gt.cells.z, Gt.cells.volumes, poro, false);
    co2_sat  = height2Sat(tmp, Gt, fluid);
    s        = [1-co2_sat, co2_sat];
    
    %% Computing hydrostatic pressure of each cell
    inplume = s(:,2) > 1/2; % All cells filled more than 50% with CO2
                            % have their cell centers inside the CO2 plume
    p = zeros(G.cells.num, 1);
    % pressure for cells in water zone is just determined by hydrostatic
    % pressure of water
    cell_z = G.cells.centroids(:,3);
    p(~inplume) = top_press + cell_z(~inplume) .* norm(gravity) * fluid.rhoWS;
    
    % pressure for cells in plume obtained by quasi-integration using cell
    % heights as steps.  
    
    % If the plume reaches the bottom, the algorithm will not work, so check
    % that first.
    assert(any(cellfun(@(x) strcmp(x, 'cartGrid'), G.type))); 
    [i,j,k] = ind2sub(G.cartDims, 1:G.cells.num);
    bottom_cells = k==max(k);
    assert(bottom_cells * double(inplume) == 0);
    
    % Do the quasi integration iteratively
    uncovered = inplume;
    zshift = prod(G.cartDims(1:2));
    while any(uncovered)
        below_ix = find(uncovered) + zshift; % indices of cells immediately below the uncovered ones
        below_targets = below_ix(~uncovered(below_ix));

        targets = below_targets - zshift;
        targets_sat = s(targets, 2) - 0.5;
        assert(all(targets_sat >= 0.0)); % otherwise, would not be in plume
        below_sat = min(s(below_targets, 2), 0.5);
        co2_frac = targets_sat + below_sat;
        dz = cell_z(below_targets) - cell_z(targets);
        assert(all(dz>0));
        pbelow = p(below_targets);
        tbelow = tinfo{1} + cell_z(below_targets) * tinfo{3} / 1000;
        dp = dz .* ((1 - co2_frac) .* fluid.rhoWS + co2_frac .* EOS.rho(pbelow, tbelow)) * norm(gravity);
        p(targets) = pbelow - dp;
        uncovered(targets) = false;
    end
end

function G = adjustVerticalCoords(G, zmin, zmax, zcellnum, exponent)
    
    tmp = linspace(0, 1, zcellnum+1);
    tmp = tmp .^ exponent;
    tmp = tmp * (zmax-zmin) + zmin;
    tmp = reshape(repmat(tmp, prod(G.cartDims(1:2)+1), 1), [],1);
    G.nodes.coords(:,3) = tmp;
end

        