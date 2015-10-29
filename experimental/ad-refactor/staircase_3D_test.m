function staircase_3D_test()

   gravity on;
   
   %% Parameters
   T = 60 + 273.15; % reservoir temperature
   p_ref = 1 * mega * Pascal;
   rhoWS = 1000;
   rhoGS = 500;
   cap_press = 0; %9.8 * 500 * 1;%0;
   duration = 1000 * year; % 100 * year
   
   %% Setting up structures
   [G, rock] = setup_grid();
   fluid = setup_fluid(p_ref, rhoWS, rhoGS, cap_press);
   initState = setupInitState(G, p_ref, rhoWS);
   model = twoPhaseGasWaterModel(G, rock, fluid, T, 0);
   schedule = initSchedule(duration);
   
   %% Running simulation
   [~, states] = simulateScheduleAD(initState, model, schedule);

   
   keyboard;
end
% ----------------------------------------------------------------------------

function [G, rock] = setup_grid()
   
   G = computeGeometry(cartGrid([3 1 2], [30 1 30]));
   rock.perm = 1 * darcy * ones(6, 1);
   rock.poro = 0.3 * ones(6,1);
end

% ----------------------------------------------------------------------------

function fluid = setup_fluid(p_ref, rhoWS, rhoGS, cap_press)

   cg = 4.6e-10;
   cw = 4.6e-10;
   
   mug = 6e-5 * Pascal * second;
   muw = 8e-4 * Pascal * second;
   mug = muw;
   
   fluid.rhoGS = rhoGS;
   fluid.rhoWS = rhoWS;
   
   fluid.rhoG = @(p, varargin) (1 + cg * (p - p_ref)) * fluid.rhoGS;
   fluid.rhoW = @(p, varargin) (1 + cw * (p - p_ref)) * fluid.rhoWS;
   
   fluid.bG = @(p, varargin) fluid.rhoG(p, varargin{:}) / fluid.rhoGS;
   fluid.bW = @(p, varargin) fluid.rhoW(p, varargin{:}) / fluid.rhoWS;
   
   fluid.muG = @(p, varargin) mug * ones(size(p)); 
   fluid.muW = @(p, varargin) muw * ones(size(p));
   
   fluid.pcGW = @(sg, p, varargin) cap_press * sg;
   
   kr = coreyPhaseRelpermAD(1, 0);
   fluid.relPerm = @(sg) deal(kr(1-sg), kr(sg));
   
end

% ----------------------------------------------------------------------------

function initState = setupInitState(G, p_ref, rhoWS)
   
   s = [1 0; 0 1; 1 0; 1 0; 1 0; 1 0];
   %s = [0.5 0.5; 0 1; 0.5 0.5; 1 0; 1 0; 1 0];
   %s = [0.3 0.7; 0 1; 0.3 0.7; 1 0; 1 0; 1 0];
   %   s = [1 0; 1 0; 1 0; 1 0; 0 1; 1 0]; % REVERSE

   initState = ...
       struct('pressure', p_ref + rhoWS * norm(gravity) * G.cells.centroids(:,3), ...
              's'       , s, ...
              'sGmax'   , s(:,2));
   
end

% ----------------------------------------------------------------------------

function schedule = initSchedule(duration)

   stepnum = 10;
   schedule = struct('control', struct('W', []) , ...
                     'step'   , struct('control', ones(stepnum, 1), ...
                                       'val', duration/stepnum * ones(stepnum, 1)));
   
end
