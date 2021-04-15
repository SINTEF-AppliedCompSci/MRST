function [myreport, mystate] = runEclipseDeck(deck)
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


  % Read deck
  mrstModule add blackoil
  mrstModule add blackoiltransport
  %deck = readEclipseDeck(deckfile_name);
  % Convert units
  deck=convertDeckUnits(deck);
  % process grid section to make grid
  G = initEclipseGrid(deck);
  G = computeGeometry(G);
  gravity on
  % make rock
  rock = initEclipseRock(deck);
  % make fluid
  fluid = initEclipseFluid(deck);
  %%-------------------------------------------------------------
  %% Make solvers at the moment this should run only for sequential splitting
  %pressure solver
  %T = computeTransECLIPSE(G, rock, deck);
  T = computeTrans(G, rock);
  psolver = @(state, w, dt, vd_tmp)               ...
            compTPFA(state, G, rock, T, fluid, dt,...
                     'wells',w,                   ...
                     'volume_correction', vd_tmp);
  clear T;
  %transport solver
  %{
  tsolver   = @(state, w,dt) ...
      explicitTransport(state,G, dt,rock, fluid,'wells',w);
  %}
  tsolver   = @(state, w,dt) ...
     implicitTransportSat(state,G, dt,rock, fluid,'wells',w);

  porvol   = poreVolume(G, rock);
  vd       = @(dt, u) ((sum(u, 2) - 1) .* porvol ./ dt);
  %--------------------------------------------------------------
  init_state = initEclipseState(G, deck, fluid)
  prev_control = nan;
  %% Define the time stepping either from ta reference run or from deck
  nnsteps=1:numel(deck.SCHEDULE.step.val);
  control = deck.SCHEDULE.step.control(1);
  W = processWells(G, rock, deck.SCHEDULE.control(control), ...
                   'InnerProduct', 'ip_tpf');
  names = {'water', 'oil', 'gas'};
  phase = isfield(deck.RUNSPEC, upper(names));
  for i=1:numel(W);
     W(i).compi=W(i).compi(phase);
  end
  init_state.wellSol=  initWellSol(W, 100*barsa);%put well pressure ??
  state={};
  t=0;
  state=init_state;
  report=[];
  mystate = {};
  t_rep  = wellCalculateProductionComp(state, W, fluid,deck.RUNSPEC, t);
  report = addToTimeStruct(report, t_rep);
  %% start simulation loop
  for step = nnsteps;
    dt = deck.SCHEDULE.step.val(step);
    if(mrstVerbose)
      disp(['Time step', num2str(convertTo(dt,day)) ]);
    end
    control = deck.SCHEDULE.step.control(step);
    if control ~= prev_control,
      W = processWells(G, rock, deck.SCHEDULE.control(control),'InnerProduct', 'ip_tpf');
      for i=1:numel(W);
         W(i).compi=W(i).compi(phase);
      end
      prev_control = control;
    end
    [c, rho, mu, u] = fluid.pvt(state.pressure, state.z); %#ok
    vd_tmp = vd(dt, u);
    state = psolver(state,W,dt,vd_tmp);
    t_rep  = wellCalculateProductionComp(state, W, fluid,deck.RUNSPEC,t);
    report= addToTimeStruct(report, t_rep);
    state = tsolver(state,W,dt);
    t = t + dt;
    mystate{end+1} = state;
    if(mrstVerbose)
      disp(['Day ', num2str(convertTo(t,day)) ]);
    end
    %plotAll(state,G);
  end
  %plotReportData(report,W)
  myreport = report;
end
%--------------------------------------------------------------------------
%-------------------------------------------------------------------------
% Plot report data
function errorrept=plotReportData(report,W)
  r_fields=fieldnames(report);
  nn=numel(r_fields);
  ncol=1;
  figure(2),clf
  for i=1:nn;
    subplot(nn,ncol,ncol*(i-1)+1),hold on
    plot(report.TIME, report.(r_fields{i}));hold on;
  end
end
%% data
function plotAll(state,G)
   figure(1),clf
   hold on
   subplot(4,1,1);hold on
   plot(1:G.cells.num, state.s(1:G.cells.num,:));
   subplot(4,1,2);hold on
   plot(1:G.cells.num, state.z(1:G.cells.num,:));
   subplot(4,1,3);hold on
   plot(1:G.cells.num, state.pressure(1:G.cells.num,:)/barsa);
   subplot(4,1,4);hold on
   plot(1:G.cells.num, state.flux(1:G.cells.num));
end


