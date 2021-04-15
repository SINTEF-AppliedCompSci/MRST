function [W, x, G, fluid, deck, report, errorrept] = ...
      runAndCompareComp(dir, only_run, error_tol)
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


  do_plot = true;

  % Read deck
  deck = readEclipseDeck(fullfile(dir, 'TEST.DATA'));

  % Convecert units
  deck = convertDeckUnits(deck);

  % process grid section to make grid

  G = processGRDECL(deck.GRID);
  G = computeGeometry(G);
  gravity reset on %off

  % make rock
  rock = initEclipseRock(deck);
  rock = compressRock(rock, G.cells.indexMap);

  % make fluid
  fluid = initEclipseFluid(deck);

  % Make solvers
  [solvers, style] = makeSolvers(G, rock, fluid, deck);
  figure(99), clf, hold on
  for kk = 1 : numel(solvers),
     plot([1 2], [kk kk], style.color{kk})
  end
  legend(style.name)
  cla, axis off


  if false,
    % initilizing from EQUIL
    init_state = initEclipseState(G, deck, fluid);
  else
    %initializing from PRESSURE and SWAT
    nph = numel(fluid.names);
    assert (nph == 2, 'Initial solution only supported for two phases.');

    sw = deck.SOLUTION.SWAT;

    init_state = struct('pressure', deck.SOLUTION.PRESSURE, ...
                        's'       , [sw, 1 - sw],           ...
                        'flux'    , zeros([G.faces.num, 1]));

    [B, B, B, B, B, B] = fluid.pvt(init_state.pressure, init_state.s); %#ok
    init_state.z = reshape(B \ reshape(init_state.s.',[],1), nph, [])';

    % init face pressures
    cellno = rldecode(1:G.cells.num, diff(G.cells.facePos), 2).';
    init_state.facePressure = ...
        accumarray(G.cells.faces(:,1), ...
                   init_state.pressure(cellno), [G.faces.num, 1]) ./ ...
        accumarray(G.cells.faces(:,1), 1, [G.faces.num, 1]);
  end

  prev_control = nan;

  %% Define the time stepping either from a reference run or from deck
  if only_run,
     % Use steping from schedule
     nnsteps = 1 : numel(deck.SCHEDULE.step.val);
  else
     % Use steping from summary eclipse output
     smry = readEclipseSummaryFmt(fullfile(dir, 'TEST'));
     tse  = smry.TIME;
     trs  = smry.TIME(smry.RptTime);

     etstep  = 2;
     nnsteps = 1 : numel(tse) - 1;
  end

  control = deck.SCHEDULE.step.control(1);
  W = processWells(G, rock, deck.SCHEDULE.control(control), ...
                   'InnerProduct', 'ip_tpf');

  t = 0;
  init_state.wellSol = initWellSol(W, 100*barsa);

  init_report = wellCalculateProductionComp(init_state, W, ...
                                            fluid, deck.RUNSPEC, t);

  state (1 : numel(solvers)) = { init_state  };  clear init_state
  report(1 : numel(solvers)) = { init_report };  clear init_report

  %% start simulation loop
  for step = nnsteps,
     if only_run,
        ddt = deck.SCHEDULE.step.val(step);
     else
        ddt = convertFrom(tse(step + 1) - tse(step), day);
     end
     dt = ddt;

     dispif(mrstVerbose, 'Time step %.2f\n', convertTo(dt, day));

     control = deck.SCHEDULE.step.control(1);
     if control ~= prev_control,
        W = processWells(G, rock, deck.SCHEDULE.control(control), ...
                         'InnerProduct', 'ip_tpf');

        for i = 1 : numel(W),
           W(i).compi = [W(i).compi(1), 1 - W(i).compi(1)];
        end

        wSol = initWellSol(W, 300*barsa);
        x.wellSol = wSol;
        prev_control = control;
     end

     for kk = 1 : numel(solvers),
        [state{kk}, report{kk}] = solvers{kk}(state{kk}, report{kk}, ...
                                              W, dt, t + dt);
     end

     t = t + dt;

     % this is done inside solvers since impes and implicite seems to have
     % different definitions ? ?
     %for kk=1:numel(solvers)
     %  t_rep  = wellCalculateProduction(state{kk}, W, fluid, t);
     %  report{kk}= addToTimeStruct(report{kk}, t_rep);
     %end

     dispif(mrstVerbose, 'Day %.2f\n', convertTo(t, day));

     if do_plot,
        if only_run,

           plotAll(state, G, style.color);

        else

           % plot and keep track of eclipse step
           etstep = plotAllEcl(dir, state, G, style.color, ...
                               etstep, t, trs, fluid);

        end
     end
  end

  plotReportData(dir, report, style.color, W, only_run)
  %err_report=checkReportData(dir,    report,W,only_run,1e-2)
end

%--------------------------------------------------------------------------

% Read Eclipse Sim
function errorrept = plotReportData(dir, report, style, W, only_run)
   if ~only_run,
      smry   = readEclipseSummaryFmt(fullfile(dir, 'TEST'));
      report = convertReportUnits2Ecl(report, smry, W);

      %find field bouth defined in report and in smry
      e_fields = fieldnames(smry);
      c_fields = e_fields(cellfun(@(n) isfield(report{1}, n), e_fields));

      ncol     = 2;
      nn       = numel(c_fields);
   else
      c_fields = fieldnames(report{1});
      nn       = numel(c_fields);
      ncol     = 1;
   end

   %NB TIME is considered a field for simplicity of code
   figure(2), clf
   for kk = 1 : numel(report),
      for i = 1 : nn,
         subplot(nn, ncol, ncol*(i-1) + 1), hold on
         plot(report{kk}.TIME, ...
              report{kk}.(c_fields{i}), [style{kk}, 'o-']); hold on
         title([c_fields{i}, ' o is mrst x is ecl'])

         if ~only_run,
            assert (ncol == 2)
            subplot(nn, ncol, ncol*(i-1) + 2), hold on

            d = report{kk}.(c_fields{i})(:, 2 : end) - ...
                smry      .(c_fields{i})(:, 2 : end);

            plot(report{kk}.TIME(2:end), d, [style{kk}, 'o-']); hold on
            title(['Diff ', c_fields{i}])

            clear d
         end
      end
   end

   % add eclipse plots in red
   if ~only_run,
      for i = 1 : nn
         subplot(nn, ncol, ncol*(i-1) + 1),hold on
         plot(smry.TIME, smry.(c_fields{i}), '-xr')
      end
   end
end

%--------------------------------------------------------------------------

function errorrept = checkReportData(dir, report, W, error_tol)
%function for checking difference from eclipse should be extended
  c_fields={'WBHP','WOPR','WWPR','WWCT'};% assumed to have right units
  nn=numel(c_fields);
  smry = readEclipseSummaryFmt(fullfile(dir, './TEST'));
  report=convertReportUnits2Ecl(report,smry,W)
  for kk=1:numel(report)
    summary_ok=true;
    for i=1:nn;
      if(isfield(smry,c_fields{i}))
        if(impes)
          nskip=10;
          disp(['Skip the ',num2str(nskip),'first points']);
          tmp_err= report.(c_fields{i})(:,(2+nskip):end-1)-smry.(c_fields{i})(:,(3+nskip):end);
          nn_norm= norm(report{kk}.(c_fields{i})(:,(2+nskip):end-1));
        else
          tmp_err= report{kk}.(c_fields{i})(:,2:end)-smry.(c_fields{i})(:,2:end);
          nn_norm= norm(report{kk}.(c_fields{i})(:,2:end));
        end
        rel_err=norm(tmp_err)/nn_norm;
        errorrept.(c_fields{i})=rel_err;
        if((rel_err>error_tol && nn_norm>0))
          summary_ok=false;
          disp(['not acceptable diff in ',c_fields{i},' value ',num2str(rel_err),' for run in dir: ',pwd])
        else
          disp(['Error in ',c_fields{i},' value ',num2str(rel_err),' for run in dir: ',pwd])
        end
      end

    end
    if(~summary_ok)
      errorrept{kk}.ok=false;
      %  error(['Error in fun at dir',pwd])
    else
      errorrept{kk}.ok=true;
    end
  end
end

%--------------------------------------------------------------------------

%% Simplest possible plotting of results
function plotAll(state,G,style)
   figure(1),clf
   hold on
   for kk=1:numel(state)
      subplot(4,1,1);hold on
      plot(1:G.cells.num, state{kk}.z(1:G.cells.num,:),[style{kk}, '.']);
      subplot(4,1,2);hold on
      plot(1:G.cells.num, state{kk}.z(1:G.cells.num,:),[style{kk}, '.']);
      subplot(4,1,3);hold on
      %plot(1:G.cells.num, state{kk}.p_ph(1:G.cells.num,:)/barsa,style{kk});
      plot(1:G.cells.num, state{kk}.pressure(1:G.cells.num,:)/barsa,style{kk});
      subplot(4,1,4);hold on
      plot(1:G.cells.num, state{kk}.flux(1:G.cells.num),style{kk});
   end
   drawnow
end

%--------------------------------------------------------------------------

function etstep = plotAllEcl(dir,   state,G,style,etstep,t,trs,fluid)
   if abs(convertTo(t, day) - trs(etstep)) < 1e-3,
      figure(1), clf, hold on

      rptfile = fullfile(dir, sprintf('TEST.F%04d', etstep - 1));
      ecl_out = readEclipseOutputFileFmt(rptfile);
      etstep  = etstep + 1;

      %{
    c_fields={'SWAT','SGAS','PRESSURE'}
    if(isfield(ecl_out,'FLOOILIp') && isfield(ecl_out,'FLOGASIp'))
      eflux=ecl_out.FLOOILIp.values+ecl_out.FLOGASIp.values;
      sat=[ecl_out.FLOOILIp.values,ecl_out.FLOGASIp.values]';
      [c, rho, mu, u, R,B] = fluid.pvt(ecl_out.PRESSURE.values*barsa, x.z);
      eflux=sum(reshape(B*sat(:),2,[]),1)';
    end
    %}

      % convert state to ecl_state for easier plotting
      % plot all state
      ecl_state = convertState(state{1}, G, fluid);
      fname = fieldnames(ecl_state);

      % find the fields that are common to both 'ecl_state' and 'ecl_out'

      p_fields = fname(cellfun(@(n) isfield(ecl_out, n), fname));

      nr   = numel(p_fields);
      ncol = 2;
      for kk = 1 : numel(state)
         ecl_state = convertState(state{kk}, G, fluid);

         for ee = 1 : numel(p_fields),
            subplot(nr, ncol, ncol*(ee-1) + 1), hold on
            v = ecl_state.(p_fields{ee})(1:G.cells.num, 1);
            plot(v, [style{kk}, '.']);

            subplot(nr, ncol, ncol*(ee-1) + 2), hold on
            d = ecl_state.(p_fields{ee}) - ecl_out.(p_fields{ee}).values;
            plot(d, [style{kk}, '.']);
         end
      end

      % add eclipse plot in red
      for ee = 1 : numel(p_fields),
         subplot(nr, ncol, ncol*(ee-1) + 1), hold on
         plot(ecl_out.(p_fields{ee}).values, 'r.')
      end

      % the flux key should be moved and implemented in convertState
      %{
    if(isfield(ecl_out,'FLOOILIp') && isfield(ecl_out,'FLOGASIp'))
      eflux=ecl_out.FLOOILIp.values+ecl_out.FLOGASIp.values;
      sat=[ecl_out.FLOOILIp.values,ecl_out.FLOGASIp.values]';
      [c, rho, mu, u, R,B] = fluid.pvt(ecl_out.PRESSURE.values*barsa, x.z);
      eflux=sum(reshape(B*sat(:),2,[]),1)';
    end
    if(isfield(ecl_out,'FLOOILIp') && isfield(ecl_out,'FLOGASIp'))
      subplot(nr,2,7),hold on
      plot(1:G.cells.num,state{kk}.flux(1:G.cells.num,1)*day);
      subplot(nr,2,8)
      plot(1:G.cells.num,x.flux(1:G.cells.num,1)*day-eflux)
    end
      %}
   end
end

%-------------------------------------------------------------------------

%% define solvers
function [solvers, style] = makeSolvers(G, rock, fluid, deck)
% define pressure solvers to test with sequential splitting

   psolvers = {};
   T = computeTransECLIPSE(G, rock, deck);
   psolvers{end+1} = @(state, w, dt, vd)     ...
      compTPFA(state, G, rock, T, fluid, dt, ... % 'bc', bc, ...
               'wells', w, ...
               'volume_correction', vd);
   clear T;

   %{
  S = computeMimeticIP(G,rock,'InnerProduct','ip_tpf');
  psolvers{end+1} = @(state,w, dt, vd)                            ...
      solveCompFlow(state, G, rock, S, fluid, dt,                ...
                    'Verbose', true,'wells',w,'volume_correction',vd);
  %'Verbose', true, 'linprob', lp,'wells',w,'volume_correction',vd);
  clear S;
  %}

   % define volume corection function
   porvol   = poreVolume(G, rock);
   vd       = @(dt, u) ((sum(u, 2) - 1) .* porvol ./ dt);

   % define transport solvers to test with sequential splitting
   tsolvers = {};
   tsolvers{end+1} = @(state, w, dt) ...
      explicitTransport(state,G, dt,rock, fluid, 'wells', w);

   %{
  tsolvers{end+1}   = @(state, w,dt) ...
      implicitTransport(state,G, dt,rock, fluid,'wells',w);
  tsolvers{end+1}   = @(state, w,dt) ...
      implicitTransportSat(state,G, dt,rock, fluid,'wells',w);
  %}

   %make all sequential splitting solvers
   solvers = cell([numel(psolvers) * numel(tsolvers), 1]);
   name    = cell([numel(psolvers) * numel(tsolvers), 1]);
   ix = 1;
   for i = 1 : numel(psolvers),
      for j = 1 : numel(tsolvers),
         solvers{ix} = @(state, report, w, dt, t)                     ...
            makeSeqSolver(psolvers{i}, tsolvers{j}, fluid, vd, state, ...
                          report, w, dt, t, deck.RUNSPEC);

         name{ix} = [func2str(psolvers{i}),' ', func2str(tsolvers{j})];

         ix = ix + 1;
      end
   end

   % Define additional solvers like impes or fully implicite to test
   %%{
   T = computeTransECLIPSE(G, rock, deck);
   solvers{end + 1} = @(state, report, w, dt, t) ...
      impesTPFA_report(state, G, T, fluid, dt, porvol, ...
                       w, report, deck.RUNSPEC, t);

   name{end + 1} = func2str(solvers{end});
   clear T;
   %}

   style.name = name;

   % give som colors for plotting
   style.color = {'b', 'k', 'g', 'c', 'y', 'm'};
end

%--------------------------------------------------------------------------

function [state, report] = ...
      impesTPFA_report(state, G, T, fluid, dt, porvol, ...
                       w, report, runspec, t)

   trans = 1 ./ accumarray(G.cells.faces(:,1), 1 ./ T, [G.faces.num, 1]);

   z0 = state.z;

   state  = impesTPFA(state, G, trans, fluid, dt, porvol, ...
                      'wells', w, 'EstimateTimeStep', false);

   z1 = state.z;
   state.z = z0;

   t_rep  = wellCalculateProductionComp(state, w, fluid, runspec, t);
   report = addToTimeStruct(report, t_rep);

   state.z = z1;
end

%--------------------------------------------------------------------------

function [state, report] = ...
      makeSeqSolver(psolver, tsolver, fluid, vd, ...
                    state, report, w, dt, t, runspec)

   [u, u, u, u] = fluid.pvt(state.pressure, state.z);  %#ok

   state  = psolver(state, w, dt, vd(dt, u));

   t_rep  = wellCalculateProductionComp(state, w, fluid, runspec, t);
   report = addToTimeStruct(report, t_rep);

   state  = tsolver(state, w, dt);
end

%--------------------------------------------------------------------------

function report = convertReportUnits2Ecl(report, smry, W)

   % define some fields and some units assume metric units for now
   u = struct('press', barsa, 'time', day, 'rate', meter^3/day, 'none', 1);

   c_fields = { 'TIME', 'WBHP' , 'WOPR', 'WGPR', 'WVPR', 'WWPR', 'WWCT' };
   c_units  = [ u.time, u.press, u.rate, u.rate, u.rate, u.rate, u.none ];

   %change order of wells to fit with eclipse
   well_id = findEclipseNumberOfWells(smry);
   ix      = cellfun(well_id, upper({ W.name })) .';

   for kk = 1 : numel(report),
      flds = find(cellfun(@(n) isfield(report{kk}, n), c_fields));

      for f = reshape(flds, 1, []),
         report{kk}.(c_fields{f}) = convertTo(report{kk}.(c_fields{f}), ...
                                              c_units(f));

         if f > 1,
            report{kk}.(c_fields{f}) = report{kk}.(c_fields{f})(ix, :);
         end

         %HACK since I know only the first two fields are not rates
         % need to be fixed using information from well structure
         if f > 2,
            report{kk}.(c_fields{f})(1,:) = abs(report{kk}.(c_fields{f})(1,:));
            report{kk}.(c_fields{f})(2,:) = report{kk}.(c_fields{f})(2,:)*0.0;
         end
      end
   end
end

%--------------------------------------------------------------------------

function ecl_state = convertState(state, G, fluid)
   ecl_state = [];

   fld_map = struct('water', 'SWAT', 'oil', 'SOIL', 'gas', 'SGAS');

   for i = 1 : numel(fluid.names),
      fld = fld_map.(lower(fluid.names{i}));
      ecl_state.(fld) = state.s(:, i);
   end

   ecl_state.PRESSURE = convertTo(state.pressure(1:G.cells.num), barsa);

   %{
% need sto be implemented
[c, rho, mu, u, R,B] = fluid.pvt(ecl_out.PRESSURE.values*barsa, x.z);
[kr] = fluid.relperm(state.s,state);
lambda=kr./mu;
%% have to use all kind of upwind trik or use the pressure
%  depends on solver FLOOILIP FLOGASIP
ecl_state.FLOVB\(x.flux(1:G.cells.num,1)
   %}
end
