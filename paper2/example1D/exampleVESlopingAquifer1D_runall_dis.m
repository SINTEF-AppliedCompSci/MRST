%% VE simulation in a standard black-oil solver
% In this example we show how to set up a standard format black-oil
% model that can be used to simulate a VE model. For the actual
% simulation, we use the fully-implicit solver in MRST from the 'ad-fi'
% module, which is based on automatic differentiation. 

try
   require deckformat ad-fi
catch% #ok<CTCH>
   mrstModule add deckformat ad-fi
end
data_dir = 'data_all_results_after/'; 
mkdir(data_dir); 
%% Parameters for the simulation
% close all

mrstVerbose true
gravity on
force_timesteps = false; 

n_fac = 10; 
T_inj = 50 * year; dt_inj = 10 * year / 5; 
T_mig = 2000 * year; dt_mig = 100 * year / 5; 

t1 = tic; 

%%
for use_dis = [true, false]; 
   for depth = [2300, 1300]
      for smooth = [false, true]
         for res_fluid = [true, false]
                                      
            [nx, ny, nz] = deal(100 * n_fac, 1, 1); % Cells in Cartsian grid
            [Lx, Ly, H]  = deal(30e3, 10e3, 50); % Physical dimensions of reservoir
            total_time   = 5 * year; % Total simulation time
            nsteps       = 40; % Number of time steps in simulation
            dt           = total_time / nsteps; % Time step length
            perm         = 1000; % Permeability in milli darcies
            phi          = 0.03; % Porosity
            
            %% Create input deck and construct grid
            % Create an input deck that can be used together with the fully-implicit
            % solver from the 'ad-fi' module. Since the grid is constructed as part of
            % setting up the input deck, we obtain it directly. 
            G = cartGrid([nx, ny, nz], [Lx, Ly, H]); 
            x = G.nodes.coords(:, 1); 
            LL = Lx * 2 / 3; 
            if(smooth)
               G.nodes.coords(:, 3) = G.nodes.coords(:, 3) + depth - LL * sin(x / LL) * tan(phi);
            else
               G.nodes.coords(:, 3) = G.nodes.coords(:, 3) + depth - LL * sin(x / LL) * tan(phi) ...
                                      + 2 * sin(2 * pi * x / 0.3e3); 
            end
            G = computeGeometry(G); 
            rock = struct('perm', 100 * milli * darcy * ones(G.cells.num, 1), ...
                          'poro', 0.2 * ones(G.cells.num, 1)); 

            %%
            W = []; 
            W = createSampleWell(W, G, rock, floor(0.1 * nx),...
                                 'Type', 'rate', 'Val', 1 * 1e6 / year,...
                                 'Radius', 0.125, 'Name', 'P1', 'Comp_i', [0 0 1]); 
            W = createSampleWell(W, G, rock, G.cartDims(1),...
                                 'Type', 'bhp', 'Val', 300 * barsa,...
                                 'Radius', 0.125, 'Name', 'P1', 'Sign',- 1, 'Comp_i', [0 1 0]); 
            
            %% Initialize data structures
            % First, we convert the input deck to SI units, which is the unit system
            % used by MRST. Second, we initialize the rock parameters from the deck; 
            % the resulting data structure may have to be post-processed to remove
            % inactive cells. Then we set up the fluid object and tell the ad-fi solver
            % that that we are working with an oil-gas system.

            % set the capillary pressure and the VE relperms explicitely
            Gt = topSurfaceGrid(G); 
            
            mu = [6e-2 * milli 8e-4] * Pascal * second; 
            rho = [760 1100] .* kilogram / meter^3; 
            if(res_fluid)
               sr = 0.21; sw = 0.11; % kwm = [0.75 0.54]; 
            else
               sr = 0; sw = 0; kwm = [1 1]; 
            end
            fluidADI = initSimpleADIFluid('mu', [mu(2) mu(2) mu(1)],...
                                          'rho', [rho(2) rho(2), rho(1)],...
                                          'n', [1 1 1]); 
            wfields = {'krO', 'krW', 'krG', 'pcOG', 'pcOW'}; 
            for i = 1:numel(wfields)
               if(isfield(fluidADI, wfields{i}))
                  fluidADI = rmfield(fluidADI, wfields{i}); 
               end
            end
            
            fluidADI.pvMultR = @(p) 1 + (1e-5 / barsa) * (p - 100 * barsa); 
            fluidADI.bO = @(p, varargin) 1 + (4.3e-5 / barsa) * (p - 100 * barsa); 
            fluidADI.BO = @(p, varargin) 1 ./ fluidADI.bO(p); 
            Temp = Gt.cells.z * 30 / 1e3 + 274 + 12; 
            p0 = Gt.cells.z(:) * norm(gravity) * 1100; 
            fluidADI.bG = boCO2(Temp, fluidADI.rhoGS); fluidADI.BG = @(p) 1 ./ fluidADI.bG(p); 

            %% 
            W(2).val = fluidADI.rhoOS * Gt.cells.z(W(2).cells) * norm(gravity); 
            % defnine relperm
            fluid = {}
            fluidADI.surface_tension = 30e-3; 
            ff_names = {'sharp interface', 'linear cap.', 'P-scaled table', 'P-K-scaled table', 'S table'}; 
            
            leg = {}; 
            rock2D = averageRock(rock, Gt); 
            opt = struct('res_gas', sr,...
                         'res_water', sw,...
                         'Gt', Gt, 'rock', rock2D); 
            results = {}; 

            for i = 1:numel(ff_names); 
               fluid = fluidADI; 
               fluid_case = ff_names{i}
               drho = (fluid.rhoOS - fluid.rhoGS); 
               switch fluid_case
                 case 'simple'
                   fluid.krG = @(sg, varargin) sg; 
                   fluid.krOG = @(so, varargin) so; 
                   fluid.pcOG = @(sg, p, varargin) norm(gravity) *       ...
                                 (fluid.rhoOS .* fluid.bO(p) -           ...
                                  fluid.rhoGS .* fluid.bG(p)) .* (sg) .* ...
                                 opt.Gt.cells.H; 
                   fluid.res_gas = 0; 
                   fluid.res_water = 0; 
                   fluid.invPc3D = @(p) 1 - (sign(p + eps) + 1) / 2; 
                   fluid.kr3D = @(s) s; 
                 case 'integrated'
                   fluid = addVERelpermIntegratedFluid(fluid, ...
                                                       'res_water'     , opt.res_water ,...
                                                       'res_gas'     , opt.res_gas ,...
                                                       'Gt'          , opt.Gt      ,...
                                                       'kr_pressure' , true        ,...
                                                       'Gt'          , opt.Gt      ,...
                                                       'int_poro'    , false       ,...
                                                       'rock'        , opt.rock); 
                   
                 case 'sharp interface'    
                   fluid = addVERelperm(fluid, Gt,...
                                        'res_water', opt.res_water,...
                                        'res_gas', opt.res_gas); 
                 case 'linear cap.'
                   fluid = addVERelpermCapLinear(fluid,...
                                                 'res_gas'   , opt.res_gas                      ,...
                                                 'res_water'   , opt.res_water                      ,...
                                                 'beta'      , 2                                ,...
                                                 'cap_scale' , 0.2 * max(opt.Gt.cells.H) *      ...
                                                               10 * (fluid.rhoOS - fluid.rhoGS) ,...  
                                                 'H'         , opt.Gt.cells.H , 'kr_pressure' , true); 
                 
                 case 'S table' 
                   C = max(opt.Gt.cells.H) * 0.4 * drho * norm(gravity); 
                   alpha = 0.5; 
                   beta = 3; 
                   samples = 100; 
                   table_co2_1d = makeVEtables('invPc3D', @(p) max((C ./ (p + C)).^(1 / alpha), opt.res_water),...
                                               'is_kscaled', false,...
                                               'kr3D', @(s) s.^beta,...
                                               'drho', drho,...
                                               'Gt', opt.Gt,...
                                               'samples', samples); 
                   S_tab = linspace(0, 1, 10)'; 
                   S_tab_w = S_tab; 
                   kr_tab_w = S_tab_w; 
                   table_water_1d = struct('S', 1 - S_tab_w, 'kr', kr_tab_w, 'h', []); 
                   fluid = addVERelperm1DTables(fluid,...
                                                'res_water', opt.res_water,...
                                                'res_gas', opt.res_gas,...
                                                'height', opt.Gt.cells.H,...
                                                'table_co2', table_co2_1d,...
                                                'table_water', table_water_1d); 
                 
                 case 'P-scaled table'
                   C = max(opt.Gt.cells.H) * 0.4 * drho * norm(gravity); 
                   alpha = 0.5; 
                   beta = 3; 
                   samples = 100; 
                   table_co2_1d = makeVEtables('invPc3D', @(p) max((C ./ (p + C)).^(1 / alpha), opt.res_water),...
                                               'is_kscaled', false, 'kr3D', @(s) s.^beta,...
                                               'drho', drho,...
                                               'Gt', opt.Gt, 'samples', samples); 
                   S_tab = linspace(0, 1, 10)'; 
                   S_tab_w = S_tab; 
                   kr_tab_w = S_tab_w; 
                   table_water_1d = struct('S', 1 - S_tab_w, 'kr', kr_tab_w, 'h', []); 
                   fluid = addVERelperm1DTablesPressure(fluid,...
                                                        'res_water', opt.res_water,...
                                                        'res_gas', opt.res_gas,...
                                                        'height', opt.Gt.cells.H,...
                                                        'table_co2', table_co2_1d,...
                                                        'table_water', table_water_1d,...
                                                        'kr_pressure', true); 
                 case 'P-K-scaled table'        
                   
                   kscale = sqrt(rock2D.poro ./ (rock2D.perm)) * fluid.surface_tension; 
                   C = 1; 
                   alpha = 0.5; 
                   beta = 3; 
                   samples = 100; 
                   table_co2_1d = makeVEtables('invPc3D', @(p) max((C ./ (p + C)).^(1 / alpha), opt.res_water),...
                                               'is_kscaled', true,....
                                               'kr3D', @(s) s.^beta,...
                                               'drho', drho,...
                                               'Gt', opt.Gt,...
                                               'samples', samples, 'kscale', kscale); 
                   S_tab = linspace(0, 1, 10)'; 
                   S_tab_w = S_tab; 
                   kr_tab_w = S_tab_w; 
                   table_water_1d = struct('S', 1 - S_tab_w, 'kr', kr_tab_w, 'h', []); 
                   fluid = addVERelperm1DTablesPressure(fluid,...
                                                        'res_water', opt.res_water,...
                                                        'res_gas', opt.res_gas,...
                                                        'height', opt.Gt.cells.H,...
                                                        'table_co2', table_co2_1d,...
                                                        'table_water', table_water_1d,...
                                                        'rock', opt.rock,...
                                                        'kr_pressure', true); 
                 otherwise
                   error('No such fluid case')
               end

               s = setupSimCompVe(Gt, rock2D); 

               if(~use_dis)
                  systemOG = initADISystemVE({'Oil', 'Gas'}, Gt, rock2D, fluid, 'simComponents', s, 'VE', true); 
               else
                  dis_rate = 5e-11; 
                  fluid.dis_rate = dis_rate; 
                  dis_max = 0.03; 
                  fluid.dis_max = dis_max; 
                  fluid.muO = @(po, rs, flag, varargin) fluidADI.muO(po); 
                  fluid.rsSat = @(po, rs, flag, varargin)   (po * 0 + 1) * dis_max; 
                              systemOG = initADISystemVE({'Oil', 'Gas', 'DisGas'}, Gt, rock2D, fluid, 'simComponents', s, 'VE', true); 
                              
               end
               
               %% Run the schedule setup in the file
               % Before we can run the schedule, we make sure that we have an initial
               % hydrostatic pressure distribution. Then we pick the schedule from the
               % input deck and start the simulator.
               % x0 = initEclipseState(G, deck, initEclipseFluid(deck)); 
               z = G.cells.centroids(:, 3); 
               clear x0; 
               % x0.pressure = ipress * barsa + (z(:)-z(end)) * norm(gravity) * fluid.rhoOS; 
               x0.pressure = W(2).val + (z(:) - z(W(2).cells)) * norm(gravity) * fluid.rhoOS; 
               x0.s(:, 1) = ones(G.cells.num, 1); 
               x0.s(:, 2) = zeros(G.cells.num, 1); 
               x0.rs = ones(G.cells.num, 1) * 0.0; 
               x0.smax = x0.s; 
               x0.smin = x0.s; 
               x0.sGmax = x0.s(:, 2); 
               dt = linspace(0.1 * year, dt_inj, 10)'; 
               dt_shift = sum(dt); 
               dt = [dt; ones(floor((T_inj - dt_shift) / dt_inj), 1) * dt_inj]; 
               dt = [dt; T_inj - sum(dt)]; 
               dt_post = linspace(0.5 * year, dt_mig, 5)'; % ' * year; 
               dt_shift = sum(dt_post); 
               dt_post = [dt_post; ones(floor((T_mig - dt_shift) / dt_mig), 1) * dt_mig]; % * year]; 
               dt_post = [dt_post; T_mig - sum(dt_post)]; 
               %% 
               control = struct('W', [], 'step', struct('val', [], 'control', [])); 
               W_post = W(2); 
               control.W = {W, W_post}; 
               control.step.val = [dt; dt_post]; 
               control.step.control = [ones(size(dt)); ones(size(dt_post)) * 2]; 
               
               systemOG.nonlinear.linesearch = false; 
               systemOG.nonlinear.maxIterations = 10; 

               systemOG.nonlinear.tol = 1e-6; 

               t2 = tic; 
               [wellSols, states] = runMrstADI(x0, Gt, systemOG, control, ...
                                               'force_step' , force_timesteps , ...
                                               'dt_min'     , 0.5 * year      , ...
                                               'report_all' , false); 
               t2 = toc(t2); 
               results{end + 1} = struct('wellSols', {wellSols}, 'states', {states}); 
            end
            zz = Gt.cells.z; 
            xc = Gt.cells.centroids(:, 1) / 1e3; 
            if(smooth)
               if(res_fluid)
                  fname = ['smooth_res_fluid', '_', '_depth_', num2str(depth)]; 
               else
                  fname = ['smooth_nores_fluid', '_', '_depth_', num2str(depth)];
               end
            else
               if(res_fluid)
                  fname = ['res_fluid','_','_depth_',num2str(depth)];
               else
                  fname = ['nores_fluid','_','_depth_',num2str(depth)];
               end
            end
            %%
            if(use_dis)
               fname = [fname,'_dis'];
            end
            if(~use_dis)
               save([data_dir,fname],'results','z','xc','control','dis_rate','dis_max','opt','t2')   
            else
               save([data_dir,fname],'results','z','xc','control','opt','t2')     
            end
         end
      end
   end
end
toc(t1)


