function [W,x,G,fluid,deck,errorrept]=runAndCompare(dir,only_run,error_tol)
%Undocumented Utility Function

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

do_plot=true;
mydir =pwd;
cd(dir);
c = onCleanup(@() cd(mydir))
% Read deck
deck = readEclipseDeck('TEST.DATA');
% Convecert units
deck=convertDeckUnits(deck);
% process grid section to make grid
G = processGRDECL(deck.GRID);
G = computeGeometry(G);
gravity on
% make rock
rock = initEclipseRock(deck);
% make fluid
fluid = initSWOFFluid('rho',deck.PROPS.DENSITY([2,1]), ...
   'mu', [deck.PROPS.PVTW(4), deck.PROPS.PVCDO(4)], ...
   'table', deck.PROPS.SWOF, 'satnum', deck.REGIONS.SATNUM );
%{
fluid = initSimpleFluid('rho',deck.PROPS.DENSITY([2,1]), ...
   'mu', [deck.PROPS.PVTW(4), deck.PROPS.PVCDO(4)],'n', [   1,   1]);
%}
% Make solvers
[solvers,style]=makeSolvers(G,rock,fluid,deck);
figure(99),clf,hold on
for kk=1:numel(solvers)
   plot([1 2],[kk kk],style.color{kk})
end
legend(style.name)
cla,axis off


if(false)
  init_state = initEclipseState(G, deck, fluid)
else
  init_state=struct('pressure',deck.SOLUTION.PRESSURE,...
                              's',deck.SOLUTION.SWAT,'flux',zeros(G.faces.num,1));
  %init_state.p_ph=[init_state.pressure(1:G.cells.num),init_state.pressure(1:G.cells.num)];
end
prev_control = nan;

%% Define the time stepping either from ta reference run or from deck
if(only_run)
   %use steping from shcedual
   nnsteps=1:numel(deck.SCHEDULE.step.val);
else
   %use steping from summary eclise file
   %[smry, rstrt] = readEclipseResults('./TEST', 'RestartData', true)
   smry = readEclipseSummaryFmt('./TEST')
   %rstrt = readEclipseRestartFmt('./TEST');
   %readEclipseResults('./TEST', 'RestartData', true);
   tse = smry.TIME;
   trs = smry.TIME(smry.RptTime);
   etstep=2;
   nnsteps=1:numel(tse)-1;
end

control = deck.SCHEDULE.step.control(1);
W = processWells(G, rock, deck.SCHEDULE.control(control), ...
   'InnerProduct', 'ip_tpf');
init_state.wellSol=  initWellSol(W, 100*barsa);
state={};
t=0;
for kk=1:numel(solvers)
  state{kk}=init_state;
  report{kk}=[];
  %[tmp,report{kk}]=solvers{kk}(state{kk},report{kk},W,0,0)
  t_rep  = wellCalculateProduction(state{kk}, W, fluid, t);
  report{kk} = addToTimeStruct(report{kk}, t_rep);
end
%% start simulation loop
for step = nnsteps;
   if(only_run)
      ddt = deck.SCHEDULE.step.val(step);
   else
      ddt = convertFrom(tse(step+1)-tse(step),day);
   end
   dt=ddt;
   if(mrstVerbose)
     disp(['Time step', num2str(convertTo(dt,day)) ]);
   end
   control = deck.SCHEDULE.step.control(1);
   if control ~= prev_control,
      W = processWells(G, rock, deck.SCHEDULE.control(control),'InnerProduct', 'ip_tpf');
      for i=1:numel(W)
         W(i).compi=[W(i).compi(1),1-W(i).compi(1)];
      end

      wSol = initWellSol(W, 300*barsa);
      x.wellSol = wSol;
      prev_control = control;
   end
   for kk=1:numel(solvers)
     [state{kk},report{kk}]=solvers{kk}(state{kk},report{kk},W,dt,t+dt);
   end
   pause(1)
   t = t + dt;
   %for kk=1:numel(solvers)
   %  t_rep  = wellCalculateProduction(state{kk}, W, fluid, t);
   %  report{kk}= addToTimeStruct(report{kk}, t_rep);
   %end

   if(mrstVerbose)
      disp(['Day ', num2str(convertTo(t,day)) ]);
   end
   if(do_plot)
      if(only_run)
         plotAll(state,G,style.color);
      else
         % plot and keep track of eclipse step
         etstep=plotAllEcl(state,G,style.color,etstep,t,trs);
      end
   end
end
%% to add the last well report (pressure is not updated where for the given
%% saturation
%for kk=1:numel(solvers)
%   impes=true;
%   if(impes)
%     [state{kk},report{kk}]=solvers{kk}(state{kk},report{kk},W,dt,t)
%   end
%end
plotReportData(report,style.color,W,only_run)
%err_report=checkReportData(report,W,only_run,1e-2)
cd(mydir)
end
%--------------------------------------------------------------------------
%-------------------------------------------------------------------------
% Read Eclipse Sim
function errorrept=plotReportData(report,style,W,only_run)
c_fields={'WBHP','WOPR','WWPR','WWCT'};% assumed to have right units
nn=numel(c_fields);
  if(~only_run)
    smry = readEclipseSummaryFmt('./TEST');
    report=convertReportUnits2Ecl(report,smry,W)
    % convert report to eclisp units
    ncol=2;
  else
     ncol=1;
  end
  figure(2),clf
  for kk=1:numel(report)
     for i=1:nn;
        subplot(nn,ncol,ncol*(i-1)+1),hold on
        plot(report{kk}.TIME, report{kk}.(c_fields{i}), ['-o',style{kk}]);hold on;
        title([c_fields{i},' o is mrst x is ecl'])
        if(~only_run)
           assert(ncol==2)
           subplot(nn,ncol,ncol*(i-1)+2),hold on
           %hack to better compeare
           %if(i>1)
           %         if(impes)
           plot(report{kk}.TIME(2:end), report{kk}.(c_fields{i})(:,2:end)-smry.(c_fields{i})(:,2:end), ['-o',style{kk}]);hold on;
           %         else
           %            plot(report{kk}.TIME(2:end), report{kk}.(c_fields{i})(:,2:end)-smry.(c_fields{i})(:,2:end), '-o');hold on;
           %         end
           title(['Diff ',c_fields{i}])
        end
     end
  end
  %% add eclpse plots
  if(~only_run)
     for i=1:nn;
        subplot(nn,ncol,ncol*(i-1)+1),hold on
        %         if(impes)
        %            plot(smry.TIME(1:end-1), smry.(c_fields{i})(:,2:end), '-x')
        %         else
        plot(smry.TIME, smry.(c_fields{i}), '-xr')
        %       end
     end
  end
end
function errorrept=checkReportData(report,W,error_tol)
 c_fields={'WBHP','WOPR','WWPR','WWCT'};% assumed to have right units
 nn=numel(c_fields);
 smry = readEclipseSummaryFmt('./TEST');
 report=convertReportUnits2Ecl(report,smry,W)
 for kk=1:numel(report)
    summary_ok=true;
    for i=1:nn;
       %plot(report.TIME, report.(c_fields{i}), '-o');hold on;
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
       %         error(['Error in fun at dir',pwd])
    else
       errorrept{kk}.ok=true;
    end
 end
end
 %% start plotting of results
function plotAll(state,G,style)
   figure(1),clf
   hold on
   for kk=1:numel(state)
      subplot(3,1,1);hold on
      plot(1:G.cells.num, state{kk}.s(1:G.cells.num,1),style{kk});
      subplot(3,1,2);hold on
      %plot(1:G.cells.num, state{kk}.p_ph(1:G.cells.num,:)/barsa,style{kk});
      plot(1:G.cells.num, state{kk}.pressure(1:G.cells.num,:)/barsa,style{kk});
      subplot(3,1,3);hold on
      plot(1:G.cells.num, state{kk}.flux(1:G.cells.num),style{kk});
   end
   drawnow
end
function etstep = plotAllEcl(state,G,style,etstep,t,trs)
if(abs(convertTo(t,day) - trs(etstep))<1e-3)
   figure(1),clf
   hold on
   %ecl_out = readEclipseOutput(sprintf('./TEST.F%04i',etstep-1));
   ecl_out = readEclipseOutputFileFmt(sprintf('./TEST.F%04i',etstep-1));
   if(isfield(ecl_out,'FLOOILI'))
      nr=3;
   else
      nr=2;
   end
   etstep=etstep+1;
   ecl_out.SWAT
   for kk=1:numel(state)
      subplot(nr,2,1),hold on
      plot(1:G.cells.num, state{kk}.s(1:G.cells.num,1),style{kk})

      plot(1:G.cells.num, state{kk}.s(1:G.cells.num,1), ...
         1:G.cells.num, ecl_out.SWAT.values(1:G.cells.num,1))
      subplot(nr,2,2),hold on
      plot(1:G.cells.num, state{kk}.s(1:G.cells.num,1)-ecl_out.SWAT.values(1:G.cells.num,1),style{kk})
      subplot(nr,2,3),hold on
      if(~isfield(ecl_out,'WAT_PRES'))
         plot(1:G.cells.num,state{kk}.pressure(1:G.cells.num,1)/barsa,style{kk});
      else
         plot(1:G.cells.num,state{kk}.p_ph(1:G.cells.num,:)/barsa,style{kk})
      end
      subplot(nr,2,4),hold on
      plot(1:G.cells.num, state{kk}.pressure(1:G.cells.num,1)/barsa-ecl_out.PRESSURE.values(1:G.cells.num,1),style{kk})
      if(isfield(ecl_out,'FLOOILI'))
         subplot(nr,2,5),hold on
         valm=state{kk}.flux(1:G.cells.num,1)*day;
         vale=ecl_out.FLOWATI.values(1:G.cells.num,1)+ecl_out.FLOOILI.values(1:G.cells.num,1);
         plot(1:G.cells.num,valm,style{kk})
         legend({'MRST','ECL'})
         subplot(nr,2,6)
         plot(1:G.cells.num,valm-vale,style{kk})
      end
      %cd(mydir)
      %return
   end
   %% add eclipse plot
   subplot(nr,2,1),hold on
   plot(1:G.cells.num,ecl_out.SWAT.values(1:G.cells.num,1),'r')
   subplot(nr,2,3),hold on
   if(~isfield(ecl_out,'WAT_PRES'))
      plot(1:G.cells.num,ecl_out.PRESSURE.values(1:G.cells.num,1),'r')
   else
      plot(1:G.cells.num,[ecl_out.PRESSURE.values(1:G.cells.num,1),ecl_out.WAT_PRES.values(1:G.cells.num,1)],'r')
   end
   if(isfield(ecl_out,'FLOOILI'))
      subplot(nr,2,5),hold on
      vale=ecl_out.FLOWATI.values(1:G.cells.num,1)+ecl_out.FLOOILI.values(1:G.cells.num,1);
      plot(1:G.cells.num,vale,'r')
   end
   drawnow
end
end
%% define solvers
function [solvers,style] = makeSolvers(G,rock,fluid,deck)
T = computeTransECLIPSE(G, rock, deck);
pc_form = 'nonwetting';
psolvers={}
%%{
psolvers{end+1} = @(state, w, dt) incompTPFA(state, G, T, fluid, 'wells',...
   w,'pc_form',pc_form,'trans_calc','upwind','well_model','sat_perf');
%}
%%{
S = computeMimeticIP(G, rock);
psolvers{end+1} = @(state, w, dt) solveIncompFlow(state, G, S, fluid, 'wells',w);
%}
%psolvers{end+1} = @(state, w, dt) ...
%   incompTPFA(state, G, T, fluid, 'wells', w,'pc_form',pc_form,'trans_calc','upwind','well_model',@wellModelPerfSat);

if(isfield(deck.RUNSPEC,'IMPES'))
   use_impl=false;
   impes=true;
else
   use_impl=true;
   impes=false;
end
%if(fullimpl)
%  full_implsolve = @(state,w,dt) iterativeFullyImplicit(state, G, T, dt, rock, fluid, 'wells', w, ...
%                                                    'Verbose', false);
%end
tsolvers={};
%%{
tsolvers{end+1} = @(state, w, dt) ...
   implicitTransport(state, G, dt, rock, fluid, 'wells', w, ...
   'Verbose', true);
%}
%%{
tsolvers{end+1} = @(state, w, dt) ...
   explicitTransport(state, G, dt, rock, fluid, 'wells', w, ...
   'Verbose', true,'ComputeDt',true);
%}
solvers={}
name={}
for i=1:numel(psolvers)
   for j=1:numel(tsolvers)
      solvers{end+1} =@(state,report,w,dt,t) makeImpesSolver(psolvers{i},tsolvers{j},fluid,state,report,w,dt,t);
      name{end+1}=[func2str(psolvers{i}),' ',func2str(tsolvers{j})];
   end
end
style.name=name;
style.color={'b','k','g','c'}
end
function [state,report] = makeImpesSolver(psolver,tsolver,fluid,state,report,w,dt,t)
  state = psolver(state,w,dt);
  t_rep  = wellCalculateProduction(state, w, fluid, t);
  report= addToTimeStruct(report, t_rep);
  state = tsolver(state,w,dt);
end
function report=convertReportUnits2Ecl(report,smry,W)
c_fields={'WBHP','WOPR','WWPR','WWCT'};%,'WVPR'}
c_units={barsa,1/day,1/day,1};%,'WVPR'}
nn=numel(c_fields);
for kk=1:numel(report)
   report{kk}.TIME=report{kk}.TIME/day;
   for i=1:nn;
      report{kk}.(c_fields{i})=report{kk}.(c_fields{i})/c_units{i};
   end
   %change order of wells to fit with eclipse
    well_id = findEclipseNumberOfWells(smry);
    ix      = cellfun(well_id, upper({ W.name })) .';
    for i=1:numel(c_fields)
      report{kk}.(c_fields{i})=report{kk}.(c_fields{i})(ix,:);
      %hack
      if(i>1)
        report{kk}.(c_fields{i})(1,:)=abs(report{kk}.(c_fields{i})(1,:));
        report{kk}.(c_fields{i})(2,:)=report{kk}.(c_fields{i})(2,:)*0.0;
      end
    end
end
end

