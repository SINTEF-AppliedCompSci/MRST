function res = setupSharpInterfaceTestCase(CO2, G, rock, varargin)

   opt = testCaseDefaults();              % load default options
   opt = merge_options(opt, varargin{:}); % override with user options
   
   % ensure any scalar parameters are expanded into vectors of the
   % appropriate size, and verify results
   num_wells = size(opt.wellPos, 1);
   opt.h0 =          expand_var(opt.h0,         prod(G.cartDims(1:2)));
   opt.wellType =    expand_var(opt.wellType,   num_wells);
   opt.schedule =    expand_var(opt.schedule,   num_wells);
   opt.bcondTypes =  expand_var(opt.bcondTypes, 4);
   
   %% Converting grid and rock to 2D
   [res.Gt, G] = topSurfaceGrid(G);
   rock2D  = averageRock(rock, res.Gt);
        
   %% To compute 's'
   res.rock = rock2D;
   
   %% Setting up fluid parameters, gravity and initial pressure
   res.gravity.theta = opt.gravityTheta;
   res.gravity.dir   = opt.gravityDir ./ norm(opt.gravityDir);
   res.fluid.water = setupWater(opt.mu(2), opt.waterDensity);
   
   % initial conditions
   res.ref_cell = opt.ref_cell; % the cell whose real depth is indep. of theta
   res.temp_grad = opt.temp_grad;
   res.top_temp = topTemp(res, opt);
   
   % initializing CO2
   ref_ix = sub2ind(res.Gt.cartDims, res.ref_cell(1), res.ref_cell(2));
   top_press = topPress(res, opt);
   ref_press = top_press(ref_ix);
   ref_temp = res.top_temp(ref_ix);
   res.fluid.CO2 = setupCO2(CO2, opt.mu(1), opt.compressible, ref_press, ref_temp);
   
   
   % (negative value for h0 signals a mass value instead, which must be
   % converted to height):
   res.initial.h0 = opt.h0;
   mc_ix = find(opt.h0 < 0); % indices of cells specified by mass rather than height
   if numel(mc_ix) > 0
       res.initial.h0(mc_ix) = columnMassToHeight(-opt.h0(mc_ix), ...
                                                  top_press(mc_ix), ...
                                                  res.top_temp(mc_ix), ...
                                                  opt.temp_grad/1000, ...
                                                  res.fluid.CO2, ...
                                                  res.fluid.water.rho, ...
                                                  res.gravity.theta, ...
                                                  res.Gt.cells.volumes(mc_ix), ...
                                                  rock2D.poro(mc_ix));
   end

   res.initial.p0 = ifacePress(res, opt);
   %[res.top_temp, res.initial.p0] = initTempAndPress(res, opt);

   
   % wells
   % @@ Following line only works for Cartesian grids.  Otherwise, check out
   % the documentaiton of 'grid_structure' and its field 'cells.indexMap'.
   warning('off', 'addWell:Sign'); % avoid annoying warning in the execution below
   num_vertical_cells = G.cartDims(3); 
   W = [];
   for i = 1:num_wells
       [px, py] = dealvec(opt.wellPos(i, :));
       % When constructing the well below, we do not set 'Val', as this can
       % change over time and is thus taken care of in the schedule
       W = verticalWell(W, G, rock, px, py, 1:num_vertical_cells, ...
                        'Type', lower(opt.wellType{i}), 'Radius', 0.1, ...
                        'Comp_i', [1 0 0]); % three components instead of two
                                            % here to avoid complaints from
                                            % initWellSolLocal (should
                                            % probably be fixed elsewhere...)
       
   end
   res.Wt = convertwellsVE(W, G, res.Gt, rock2D);
   
   warning('on', 'addWell:Sign'); % avoid annoying warning in the execution below
   
   %% Boundary conditions (needed by 'system')
   res.bc = [];
   orient = {'LEFT', 'RIGHT', 'BACK', 'FRONT'};
   for i = 1:4
       btype = upper(opt.bcondTypes{i});
       switch btype
         case 'PRESSURE'
           rhoW = res.fluid.water.rho;
           fpr = facePressures(res.Gt, res.initial.p0, orient{i}, res.gravity, rhoW);
           res.bc = pside(res.bc, res.Gt, orient{i}, fpr);
                          
         case 'FLUX'
           % @@ Currently, 'FLUX' is associated with 'no-flow'.  This can be
           % generalized later.
           res.bc = fluxside(res.bc, res.Gt, orient{i}, 0);
         otherwise 
           error('Unknown boundary condition type specified');
       end
   end
   
   %% looping information

   res.numSteps = opt.timesteps;
   res.dt = opt.simTime/opt.timesteps;
   res.schedule = setupSchedule(opt.schedule);

end

%%                LOCAL HELPER FUNCTIONS                               
%%                                                                     

function depth = computeRealDepth(Gt, gravity, ref_cell)
    % compute the depth, taking the inclination of the aquifer into account

    ref_ix = sub2ind(Gt.cartDims, ref_cell(1), ref_cell(2));
    shift = bsxfun(@minus, Gt.cells.centroids, Gt.cells.centroids(ref_ix, :));
    shift = shift * gravity.dir(:);
    depth = Gt.cells.z(ref_ix) + ...
            (Gt.cells.z - Gt.cells.z(ref_ix)) * cos(gravity.theta) - ...
            shift * sin(gravity.theta);
    
end

%%                                                                              
function temp = topTemp(res, opt)
    depth_diff = computeRealDepth(res.Gt, res.gravity, res.ref_cell) - opt.ref_depth;
    temp = opt.ref_temp + (depth_diff .* res.temp_grad ./ 1000);
end

%%                                                                              
function press = topPress(res, opt)
    depth_diff = computeRealDepth(res.Gt, res.gravity, res.ref_cell) - opt.ref_depth;
    press = opt.ref_press + norm(gravity) * res.fluid.water.rho * depth_diff;
end

%%                                                                              
function press = ifacePress(res, opt)
    g_cos_t = norm(gravity) * cos(res.gravity.theta);
    press = topPress(res, opt) + g_cos_t * res.fluid.water.rho * res.initial.h0;
end

% %%                                                                              
% function [top_temp, iface_press] = initTempAndPress(res, opt)
%     real_depth = computeRealDepth(res.Gt, res.gravity, res.ref_cell);
%     depth_diff = real_depth - opt.ref_depth;

%     top_temp = opt.ref_temp + (depth_diff .* res.temp_grad ./ 1000);
%     iface_press = opt.ref_press + ...
%                   norm(gravity) * res.fluid.water.rho * (depth_diff + res.initial.h0 * cos(res.gravity.theta));
% end
%%                                                                              

function res = setupSchedule(schedule)
% Produces a schedule on the form of a 2D table where rows correspond to
% events and columns represent: timestep | well_ix | value
    num_wells = numel(schedule);
    res = [];
    for i = 1:num_wells
        well_plan = schedule{i}; % matrix with 2 columns (timestep and val)
        num_events = size(well_plan, 1);
        well_plan = [well_plan(:,1), i * ones(num_events,1), well_plan(:,2)];
        res = [res; well_plan];
    end
    % sort the resulting schedule according to timestep
    res = sortrows(res, 1);
end
%%                                                                              
function pface = facePressures(Gt, p, side, grav, rhoW)
           
    ix = boundaryFaceIndices(Gt, side, [1:Gt.cartDims(1)], [1:Gt.cartDims(2)], []);
    % identifying indices of the cells containing these faces
    assert(prod(double(Gt.faces.neighbors(ix,:)), 2) == 0); % 
    cell_ix = sum(Gt.faces.neighbors(ix, :), 2);

    dh = -(Gt.faces.centroids(ix,:) - Gt.cells.centroids(cell_ix,:)) * grav.dir' ...
         * sin(grav.theta);
    dp = rhoW * norm(gravity) * cos(grav.theta) * dh;
        
    pface = p(cell_ix) + dp;
end
%%                                                                              

function water = setupWater(mu, rho)
% This function can be elaborated later, if water properties are to be
% governed by a proper EOS.  For now, water properties are constant.
    water.mu = mu;
    water.rho = rho;
end
%%                                                                              
function CO2fluid = setupCO2(CO2Density, mu, compressible, ref_press, ref_temp)
    CO2fluid.mu = mu;
    CO2fluid.compressible = upper(compressible);
    
    if (strcmpi(compressible, 'INCOMPRESSIBLE'))
        const_rho = CO2Density.rho(ref_press, ref_temp);
        CO2fluid.rho    = @(p,t) ones (num_elements(p), 1) * const_rho;
        CO2fluid.beta   = @(p,t) zeros(num_elements(p), 1);
        CO2fluid.gamma  = @(p,t) zeros(num_elements(p), 1);
        CO2fluid.beta2  = @(p,t) zeros(num_elements(p), 1);
        CO2fluid.gamma2 = @(p,t) zeros(num_elements(p), 1);
        CO2fluid.chi    = @(p,t) zeros(num_elements(p), 1);
    else
        CO2fluid.rho    = @CO2Density.rho;
        CO2fluid.beta   = @CO2Density.beta;     
        CO2fluid.gamma  = @CO2Density.gamma;        
        CO2fluid.beta2  = @(p,t) CO2Density.rhoDPP(p,t) ./ CO2Density.rho(p,t);
        CO2fluid.gamma2 = @(p,t) CO2Density.rhoDTT(p,t) ./ CO2Density.rho(p,t);        
        CO2fluid.chi    = @(p,t) CO2Density.rhoDPT(p,t) ./ CO2Density.rho(p,t);
    end
        
end
%%                                                                              
function res = num_elements(x)
    
    if isa(x, 'ADI')
        res = numel(x.val);
    else
        res = numel(x);
    end    
end

%%                                                                               
function h = columnMassToHeight(masses, Pcap, Tcap, Tgrad, CO2rhofun, rhoWater, theta, col_area, poro)
    N = numel(masses); 
    for i = 1:N
        h(i) = columnHeight(Pcap(i), Tcap(i), Tgrad, CO2rhofun, rhoWater, theta, masses(i), col_area(i), poro(i), true);
    end
end
