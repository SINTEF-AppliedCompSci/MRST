function [myreport, mystate] = runEclipseDeckIncomp(deck, varargin)
%{
Copyright 2009-2021 SINTEF Digital, Mathematics & Cybernetics.

This file is part of The MATLAB Reservoir Simulation Toolbox (MRST).

MRST is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

MRST is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with MRST.  If not, see <http://www.gnu.org/licenses/>.
%}


   opt = struct('use_mimetic', false, 'use_reorder', true);
   opt = merge_options(opt, varargin{:});

  % Read deck
  % Convert units
  deck=convertDeckUnits(deck);
  % process grid section to make grid
  G=initEclipseGrid(deck);
  G=mcomputeGeometry(G);
  gravity on
  % make rock
  rock = initEclipseRock(deck);
  % make fluid
  fluid_ecl = initEclipseFluid(deck);
  mu        = [deck.PROPS.PVTW(1,4), deck.PROPS.PVCDO(1,4)];
  density   = deck.PROPS.DENSITY(1:2);
  alpha     = [2, 2];
  fluid     = initSimpleFluid('mu' , mu     , ...
                              'rho', density, ...
                              'n'  , alpha);

  %%-------------------------------------------------------------
  %% Make solvers at the moment this should run only for sequential splitting
  %pressure solver
  %%
  %T = computeTransECLIPSE(G, rock, deck);
  T = computeTrans(G, rock);
  linsolve_p = @mldivide; %@(S, h) agmg(S, h,  1, 5.0e-11, 1000, 0);
  linsolve_t = @mldivide; %@(J, F) agmg(J, F, 50, 5.0e-11, 2000, 0);

  if opt.use_mimetic,
     psolver = @(x,W) ...
        solveIncompFlow(x, G, S, fluid, 'wells', W, 'LinSolve', linsolve_p);
  else
     psolver = @(x,W) ...
      incompTPFA(x, G, T, fluid, 'wells', W, 'LinSolve', linsolve_p);
  end

  if ~opt.use_reorder
     tsolver = @(x, W, dt) ...
        implicitTransport(x, G, dt, rock, fluid, 'wells', W, ...
        'LinSolve', linsolve_t);
  else
     require mex/reorder-C

     fluid.param = struct('viscw' , mu(1)   , 'visco', mu(2)   , ...
                          'srw'   , 0       , 'sro'  , 0       , ...
                          'nw'    , alpha(1), 'no'   , alpha(2), ...
                   ...
                          'satnum', ones([G.cells.num, 1], 'int32'));

     tsolver = @(x, W, dt) ...
        implicitTransportReorder(x, G, dt,rock, fluid, 'wells', W);
  end

  %--------------------------------------------------------------
  init_state = initEclipseState(G, deck, fluid_ecl);
  prev_control = nan;
  %% Define the time stepping either from ta reference run or from deck
  nnsteps=1:numel(deck.SCHEDULE.step.val);
  control = deck.SCHEDULE.step.control(1);
  W = processWells(G, rock, deck.SCHEDULE.control(control), ...
                   'InnerProduct', 'ip_tpf');
  init_state.wellSol=  initWellSol(W, 100*barsa);%put well pressure ??
  t=0;
  state=init_state;
  t_rep  = wellCalculateProduction(state, W,fluid, t);
  report=[];
  mystate = {};
  report = addToTimeStruct(report, t_rep);
  %% start simulation loop
  p_t=0;t_t=0;

  ndigits = floor(log10(nnsteps(end))) + 1;

  for step = nnsteps;
    dt = deck.SCHEDULE.step.val(step);

    dispif(mrstVerbose, 'Time step = %.3e\n', convertTo(dt, day));

    control = deck.SCHEDULE.step.control(step);
    if control ~= prev_control,
      W = processWells(G, rock, deck.SCHEDULE.control(control), ...
                       'InnerProduct', 'ip_tpf');

      % Use first two components of injection composition only.
      compi     = arrayfun(@(w) w.compi(1:2), W, 'UniformOutput', false);
      [W.compi] = compi{:};

      % RESV wells are 'rate' wells in MRST parlance.
      [W(strcmp({W.type}, 'resv')).type] = deal('rate');

      prev_control = control;
    end
    t0 = tic; state = psolver(state,W);  dt_cpu = toc(t0);
    fprintf('[%0*d]: Pressure:  %12.5f [s]\n', ndigits, step, dt_cpu);
    p_t=p_t+dt_cpu;

    t_rep  = wellCalculateProduction(state, W, fluid,t);
    report= addToTimeStruct(report, t_rep);
    t0 = tic; state = tsolver(state,W,dt);dt_cpu = toc(t0);
    fprintf('[%0*d]: Transport: %12.5f [s]\n', ndigits, step, dt_cpu);
    t_t=t_t+dt_cpu;

    t = t + dt;
    mystate{end+1} = state;
    dispif(mrstVerbose, 'Day = %.3e\n', convertTo(t, day));
  end

  blank = repmat(' ', [ndigits - 1, 1]);
  fprintf('Sum%s: Pressure:  %12.5f [s]\n', blank, p_t);
  fprintf('Sum%s: Transport: %12.5f [s]\n', blank, t_t);

  myreport = report;
end
