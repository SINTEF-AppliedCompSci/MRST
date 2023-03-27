function [Gt, optim] = verySimpleExample(varargin)

   opt.res = 33;
   opt.timesteps = 4;
   opt.min_rate = 0.5;
   opt.max_rate = 10;
   opt.timestep_len = 1 * day;
   opt = merge_options(opt, varargin{:});
   
   moduleCheck('ad-core');
   gravity on;
   
   ref_pressure = 100 * barsa;
   rate = 1; %0.18;%5;%sqrt(eps);
   
   G = computeGeometry(cartGrid([opt.res opt.res 1], [1 1 0.02] * kilo * meter));
   Gt = topSurfaceGrid(G);
   rock.perm = 100 * milli * darcy * ones(Gt.cells.num, 1);
   rock.poro = 0.2 * ones(Gt.cells.num, 1);
   rock = averageRock(rock, Gt);
   
   fluid = makeVEFluid(Gt, rock, 'simple', ...
                       'fixedT', 40 + 275.25, ...
                       'co2_rho_pvt', [], ...
                       'wat_rho_pvt', [], ...
                       'pvMult_fac',  0);
                       
   wcell = (opt.res^2 + 1)/2; %5;
   W = addWell([], Gt, rock, wcell    , ...
               'Type'   , 'rate'      , ...
               'Val'    , rate        , ...
               'Radius' , 0.3         , ...
               'Comp_i' , [0, 1]      , ...
               'name'   , ['I']);
   % W = addWell(W, Gt, rock, 1    , ...
   %             'Type'   , 'rate'      , ...
   %             'Val'    , rate        , ...
   %             'Radius' , 0.3         , ...
   %             'Comp_i' , [0, 1]      , ...
   %             'name'   , ['I']);
   
   schedule.control.W = W;
   bfaces = find(any(Gt.faces.neighbors==0,2)); % identify boundary faces
   schedule.control.bc = addBC([], bfaces, 'pressure', ref_pressure, ...
                                   'sat', [1 0]);
   schedule.step.val = ones(1, opt.timesteps) * opt.timestep_len;
   schedule.step.control = ones(1, opt.timesteps);
                                   
   
   schedule.control(2).W = W;
   schedule.control(2).W.val = opt.min_rate;
   schedule.control(2).bc = schedule.control(1).bc;
   schedule.step.control(2:end) = 2;
   
   initState.pressure = ref_pressure * ones(Gt.cells.num, 1);
   initState.s = repmat([1 0], Gt.cells.num, 1);
   initState.sGmax = initState.s(:,2);
   
   model = CO2VEBlackOilTypeModel(Gt, rock, fluid);
   min_rates = opt.min_rate * ones(numel(W), 1);
   max_rates = opt.max_rate   * ones(numel(W), 1);
   
   
   [optim, init, history] = ...
       optimizeRates(initState, model, schedule, min_rates, max_rates, ...
                     'last_control_is_migration', false);
   
end