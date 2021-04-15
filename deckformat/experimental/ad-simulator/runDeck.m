function varargout = runDeck(deckfile,varargin)
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

%mrstModule add ad-blackoil ad-core mrst-gui ad-props deckformat
%try
%linearsolver=struct('eliptic_solver','agmg','method','CPR');
%linearsolver=struct('eliptic_solver','backslash','method','backslash');

linearsolver=struct('method','CPR_agmcl');
opt=struct('method','from_deck',...
    'linearsolver',linearsolver,...
    'test_run',false,...
    'resultdir',[],...
    'outputdir','output',...
    'force_timestep',false,...
    'nonlinear','olympus',...
    'schedule_from_output',false);
opt=merge_options(opt,varargin{:});
opt.want_results = nargout>0;
%% Set up model
%%
[dirname,filename, ext] = fileparts(deckfile);
%name=deckfile(1:end-5);
outputdir=opt.outputdir;
%a=split(name,'/');
%dirname=char(a(end));
deck = readEclipseDeck(deckfile);
deck = convertDeckUnits(deck);

gravity reset on
switch opt.method
    case 'from_deck'
        
        G = initEclipseGrid(deck);
        G = computeGeometry(G);
        
        rock  = initEclipseRock(deck);
        rock  = compressRock(rock, G.cells.indexMap);
        
        % Create a special ADI fluid which can produce differentiated fluid
        % properties.
        fluid = initDeckADIFluid(deck);
        model = selectModelFromDeck(G, rock, fluid, deck);
        %%
        if(isfield(deck.SOLUTION,'EQUIL'))
            regions = getInitializationRegionsDeck(model, deck);
            [state0, pressures] = initStateBlackOilAD(model, regions);
        else
            state0 = solutionToMrstState(G,deck.SOLUTION)
        end
        %%
        model = selectModelFromDeck(G, rock, fluid, deck);
    case 'init_output'
        % do not work for opm since RSPEC is not written
        if(isempty(opt.resultdir))
            casenm = fullfile(dirname,filename);
        else
            casenm  = fullfile(opt.resultdir,filename);
        end
        init = readEclipseOutputFileUnFmt([casenm,'.INIT']);
        grid = readEclipseOutputFileUnFmt([casenm,'.EGRID']);
        [fluid, pvtDeck] = initFluidFromOutput(init);
        %[G, rock, N, T] = eclOut2mrst(init, grid); use sim-grid beacuse
        % non-matching mrst-grid du to pinch gridpinch
        [G, Gs, rock, info] = eclOutToGrids(init, grid, 'outputRock', false, 'outputSimGrid', true);
        T = info.T;
        %[G, rock, N, T] = eclOut2mrst(init, grid);
        warning('May not work if MRST grid and OPM grid is different i,e aquifers +++');
        %states = convertRestartToStates(casenm, G,'steps',1,'use_opm',true);
        [states, rstrt] = convertRestartToStates(casenm, Gs);
        %wellSols= convertReservoirFluxesToSurface(casenm,G,'steps',1)
        
        % set welsols
        %error('Have to put correctwells wellSols')
        %fluid = initDeckADIFluid(deck);
        %model = selectModelFromDeck(Gs, [], 'trans', T, 'porv', info.PORV);
        model = selectModelFromDeck(Gs, [], fluid, deck);
        model.operators = setupOperatorsTPFA(Gs, rock, 'trans', T, 'porv', info.PORV);
        %model.operators=setupOperatorsTPFA(G, [], 'neighbors', N , 'trans', T, 'porv', G.PORV);
        model.fluid = fluid;
        dispif(mrstVerbose, 'Computing and adding fluxes ...\n');
        %states = addFluxOilWater(states, model);
        state0=states{1};
        state0 = rmfield(state0,'wellSol')
        %model.setupOperators   
        %error('To be implemented')
        G=Gs;
    otherwise
        error('Valid methods are: from_deck,')
end
ncomp = model.water + model.oil + model.gas;

%schedule = convertDeckScheduleToMRST(model, deck);
% model.AutoDiffBackend = SparseAutoDiffBackend();
model.AutoDiffBackend = DiagonalAutoDiffBackend();
model.FacilityModel = UniformFacilityModel(model);
%%
% Set maximum limits on saturation, Rs and pressure changes
%model.drsMaxRel = .2;
%model.dpMaxRel  = .2;
%model.dsMaxAbs  = .05;

if(~opt.schedule_from_output)
    schedule = convertDeckScheduleToMRST(model, deck);
    schedule_org=schedule;
else
    %deck = convertDeckUnits( dec );
    %schedule = convertDeckScheduleToMRST(model, deck.SCHEDULE);
    % set directly
    schedule = deck.SCHEDULE;
    assert(all(schedule.step.control==1));
    if(isempty(opt.resultdir))
            casenm = fullfile(dirname,filename);
    else
            casenm  = fullfile(opt.resultdir,filename);
     end
    [states, rstrt] = convertRestartToStates(casenm, G);
    W = states{end}.wellSol; % don't set to first state! 
    % remove "third phase"
    for k = 1:numel(W)
        W(k).compi = W(k).compi(1:2);
        W(k).lims  = [];
    end
    schedule.control = struct('W', W);  
    % use same ministeps as ECLIPSE
    smry = readEclipseSummaryUnFmt(casenm);
    dt = diff(smry.get(':+:+:+:+', 'TIME', ':'));
    schedule.step.val = dt(:)*day;
    schedule.step.control = ones(numel(dt), 1); 
end
% sett up linear solver

switch opt.linearsolver.method
    case 'CPR'
        if model.gas && (isprop(model, 'disgas') && model.disgas) && ...
                (isprop(model, 'vapoil') && model.vapoil)
            
            
            linsolve = CPRSolverAD('ellipticSolver', pressureSolver);
        else
            linsolve = SimpleCPRAD('ellipticSolver', pressureSolver);
        end
        switch opt.linearsolver.eliptic_solver
            case 'agmg'
            try
                mrstModule add agmg
            pressureSolver = AGMGSolverAD('tolerance', 1e-4);
        catch
            pressureSolver = BackslashSolverAD();
        end
    case 'amgcl'
        pressureSolver = AMGCLSolverAD('coarsening', 'smoothed-aggregation',...
                    'relaxation', 'spai0', ...
                    'solver', 'bicgstab', ...
                    'preconditioner', 'amg')
        
    case 'backslash'
        pressureSolver = BackslashSolverAD();
    otherwise
        error('valid solver types: agma backslash')
        end
        
        
    case 'CPR_agmcl'
        linsolve = AMGCL_CPRSolverAD('block_size', size(state0.s,2), 'maxIterations', 150, 'tolerance', 1e-5);
        % Optimal?
        linsolve.setCoarsening('aggregation')
        linsolve.setRelaxation('ilu0');        
        linsolve.setSRelaxation('ilu0');
        linsolve.amgcl_setup.verbose = false;      
        linsolve.amgcl_setup.npre = 1;
        linsolve.amgcl_setup.npost = 1;
        linsolve.amgcl_setup.ncycle = 1;        
        linsolve.amgcl_setup.direct_coarse = false;
        linsolve.amgcl_setup.drs_eps_dd = 0.3;
        linsolve.amgcl_setup.drs_eps_ps = 0.1*linsolve.amgcl_setup.drs_eps_dd;
        linsolve.doApplyScalingCPR = true;
        linsolve.trueIMPES = false;
        linsolve.amgcl_setup.use_drs = true;
        %linsolve.amgcl_setup.use_drs = true;
        linsolve.verbose = true;
        if isa(model.AutoDiffBackend, 'DiagonalAutoDiffBackend')
            %lsolve.reduceToCell = true;
            linsolve.keepNumber = G.cells.num*ncomp;
            linsolve.amgcl_setup.active_rows = G.cells.num*ncomp;
        end
        linsolve.verbose = true;
        linsolve.amgcl_setup.verbose = true;
        
    case 'backslash'
        linsolve = BackslashSolverAD();
    case 'amgcl_ilu'
        lsolve = AMGCLSolverAD('preconditioner', 'relaxation', 'relaxation', 'ilu0', 'maxIterations', 150, 'solver', 'idrs');
    %     lsolve.applyRightDiagonalScaling = false;
        lsolve.tolerance = 5e-3;
        ndof = ncomp*G.cells.num;
        ordering = getCellMajorReordering(G.cells.num, ncomp, ndof);
        lsolve.variableOrdering = ordering;
        lsolve.equationOrdering = ordering;
        
    case 'amgcl'      
        
        linsolve = AMGCLSolverAD('coarsening', 'aggregation',...
                    'relaxation', 'ilut', ...
                    'solver', 'gmres', ...
                    'preconditioner', 'relaxation');
        
    otherwise
        error('Valid solvers is backslash')
end
if isa(model.AutoDiffBackend, 'DiagonalAutoDiffBackend')
    linsolve.reduceToCell = true;
    linsolve.keepNumber = G.cells.num*ncomp;
    linsolve.amgcl_setup.active_rows = G.cells.num*ncomp;
end
%linsolve = BackslashSolverAD()

%%
nstep=3;
if(opt.test_run)
    schedule=schedule_org;
    schedule.step.control = schedule_org.step.control(1:nstep)
    schedule.step.val = schedule_org.step.val(1:nstep)
    mrstVerbose true;
end
%%
% define output
%fn = getPlotAfterStep(state0, model, schedule, ...
%    'plotWell', false, 'plotReservoir', false);

%if(~opt.want_results)
  [rph, oph] = getRunDeckOutPUT(outputdir, filename);
%else
%  rph=[];
%  oph=[];
%end

reportfile=fullfile(outputdir,filename,'initstate');
mkdir(fullfile(outputdir,filename))
save(reportfile,'state0');
simtime = tic();
failed=false;
%try
%%
%tsel = StateChangeTimeStepSelector('targetProps', {'s', 'pressure'}, 'targetChangeAbs', [0.3, 40*barsa], ...
%    'targetChangeRel', [inf, 0.5], 'firstRampupStepRelative', 0.1);


switch opt.nonlinear
    case 'normal'
    nls=NonLinearSolver();
    if(~opt.force_timestep)
        tsel = StateChangeTimeStepSelector('targetProps', {'s', 'pressure'}, 'targetChangeAbs', [0.3, 40*barsa], ...
        'targetChangeRel', [inf, 0.5], 'firstRampupStepRelative', 0.1)
        nls.timeStepSelector = tsel;
    end
    nls.LinearSolver = linsolve;
    case 'olympus'
        nls = NonLinearSolver();
        nls.timeStepSelector = IterationCountTimeStepSelector('firstRampupStep',... 
                    1*day, 'targetIterationCount', 5,...
                    'maxTimestep', 30*day);
        nls.maxIterations = 12;
        nls.LinearSolver = linsolve;
        nls.useRelaxation = true;
        model.toleranceCNV = 1e-2;
        model.toleranceMB = 1e-5;
        %model.toleranceCNV = 1e-3;
        %model.toleranceMB = 1e-7q
        ;
    otherwise
        error('nonlinear not set');
end


if(opt.want_results)
    [wellsols, states, reports] =...
        simulateScheduleAD(state0, model, schedule, ...
        'NonLinearSolver', nls,...
        'OutputHandler',oph,...
        'ReportHandler', rph,...
        'afterStepFn',     [], ...
        'restartStep',     1,...
        'OutputMinisteps', true);
    results={wellsols,states,reports ,model}
    for i = 1:nargout
        varargout{i} = results{i};
    end
else
    % to keep maximum numer of ministeps in memory
    simulateScheduleAD(state0, model, schedule, ...
        'NonLinearSolver', nls,...
        'OutputHandler',oph,...
        'ReportHandler', rph,...
        'afterStepFn',     [], ...
        'restartStep',     1,...
        'OutputMinisteps', true);
end
%catch
%failed=true;
reports={};
%end
simtime=toc(simtime);
% save report
reportfile=fullfile(outputdir,filename,'report');
save(reportfile,'reports','failed','simtime','opt','deckfile','schedule');
%catch
%    disp('Running failed')
%end
end
function [rph, oph] = getRunDeckOutPUT(outputdir, dirname)
    rph= ResultHandler('writeToDisk',true,...
                       'storeInMemory',false,...
                       'dataDirectory',outputdir,...
                       'dataPrefix','report',...
                       'dataFolder',dirname,...
                       'saveflags','',...
                       'cleardir',false,...
                       'verbose',true);

    oph= ResultHandler('writeToDisk',true,...
                       'storeInMemory',false,...
                       'dataDirectory',outputdir,...
                       'dataPrefix','states',...
                       'dataFolder',dirname,...
                       'saveflags','',...
                       'cleardir',false,...
                       'verbose',true);
end
