function model = Run_timesteptest(model)
    static = SetupStaticParams(model); % set static struct
    model.static = static;
    state0 = Initialize(model); % set initial condition
    state0.facePressure = zeros(model.grid.G.faces.num,1) * model.experiment.schedule.pout.value;
    model.state0 = state0;
    model.state = state0;
    dynamic = SetupDynamicParams(model); % initialize dynamic struct 
    model.dynamic = dynamic;
    rampupsteps = model.simulation.rampupsteps.value;
    %---------------------------------------  
    scheduleTable = model.experiment.schedule.procedure;
    dt = model.static.dt;
    %---------------------------------------
    G     = model.grid.G;
    rock  = model.rock;
    fluid = model.fluid;
    model.twoPhaseOilWaterModel = TwoPhaseOilWaterModel(G, rock, fluid);
    if not(strcmp(model.simulation.type,strcat('Historymatch')))
        fprintf('Simulation of the schedule is started.\n')
    end
    %---------------------------------------
    for i = 1 : height(scheduleTable)
        params = dynamic.params;
        params.periodStart = [params.periodStart,scheduleTable{i,1}];
        params.periodEnd = [params.periodEnd,scheduleTable{i,2}];
        params.periodInterval = params.periodEnd - params.periodStart;
        scheduleRow = scheduleTable{i,:};
        model.bc = BoudaryConditions(model,scheduleRow);
        processType = model.experiment.process.type; 
        if(strcmpi(processType,'cent'))
            model.gravity = scheduleRow(3);
            model = SetupCentrifuge(model);
        else
            gravity off
        end
        twoPhaseOilWaterModel = model.twoPhaseOilWaterModel;
        %---------------------------------------
        periodInterval = params.periodInterval(end);
        n = periodInterval/dt;
        rampup = rampupTimesteps(periodInterval, periodInterval/n,rampupsteps);
        schedule = simpleSchedule(rampup, 'bc', model.bc);
        schedule_small = compressSchedule(schedule);
        
%         seqModel = getSequentialModelFromFI(twoPhaseOilWaterModel);
%         stepSel = StateChangeTimeStepSelector(...
%           'targetProps', 'pressure',...            % Saturation as change target
%           'targetChangeAbs', 100,...       % Target change of 0.01
%           'firstRampupStepRelative', 0.01); % Initial rampup step is dt0/100
%         seqModel.transportNonLinearSolver.timeStepSelector = stepSel;
        
        rampup = 1*minute;
        targetIts = 20;
        timestepper = ...
           IterationCountTimeStepSelector('targetIterationCount', targetIts,...
                                          'minRelativeAdjustment', sqrt(eps),...
                                          'maxRelativeAdjustment', 100, ...
                                          'firstRampupStep',       rampup, ...
                                          'verbose', true);

        % Instantiate a nonlinear solver with the timestep class as a
        % construction argument.
        nonlinear = NonLinearSolver('timeStepSelector', timestepper, ...
                                    'maxiterations', targetIts);

        % simulate the schedule
        %---------------------------------------
        verbose = true;
        schedule.counter = i;
        [~, model.state, model.schedulereport] = simulateScheduleAD(...
                    model.state, twoPhaseOilWaterModel, schedule_small, ...
                    'Verbose', verbose, ...
                     'nonlinearSolver', nonlinear, 'outputMinisteps', true);
%         [model.state, model.schedulereport] = Simulate(...
%                       model.state,twoPhaseOilWaterModel,schedule_small,...
%                       verbose);
%         [~, model.state, model.schedulereport] = evalc...
%             ('Simulate(model.state,twoPhaseOilWaterModel,schedule, verbose)');
        %---------------------------------------
        params.scheduleSteps = rampup;
        params.cumScheduleSteps = [params.cumScheduleSteps;...
        cumsum(params.scheduleSteps) + params.cumScheduleSteps(end)];
        params.schedule = [params.schedule; schedule];
        params.fromIdx = params.toIdx;
        params.counter = params.counter + length(schedule.step.control);        
        params.toIdx   = params.counter;
        params.plot = false;
        dynamic.params = params;             
        dynamic.states = [dynamic.states;model.state]; 
        dynamic.schedulereport = [dynamic.schedulereport; model.schedulereport];
        model.dynamic  = dynamic;
        dynamic = UpdateDynamicParams(model); % update dynamic parameters
        model.dynamic = dynamic;
        if not(strcmpi(model.simulation.type,strcat('Historymatch')))
            params.plotNo = params.plotNo + 1; 
            Visualize(model); 
            frames(params.plotNo) = getframe(gcf); 
        end
        dynamic.params.plotNo = params.plotNo;
    end    
    if not(strcmp(model.simulation.type,strcat('Historymatch')))
        model.dynamic.params.frames = frames;
%         if strcmp(model.plot.style,'normal')
%             f = findobj('type','figure','name','Visualization');
%             f.WindowState = 'minimized';
%         end
        if not(isempty(model.output))
            if(model.output.include)
                SaveResults(model);
            end
        end
    end
%---------------------------------------
% saving parameters for history match
%---------------------------------------
    model = save_for_HM(model,state0);
end