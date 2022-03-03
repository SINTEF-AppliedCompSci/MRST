function model = Run(model, visualization)
% <keywords>
%
% Purpose : run the simulation schedule
%
% Syntax :
%   model = Run(model, visualization)
%
% Input Parameters :
%   model: struct containing the modeling specifications
%   visualization: boolean controling the visualization during the
%   simulation
%
% Return Parameters :
%   model: struct containing the simulation results
%
% Description :
%
% Author : 
%    Siroos Azizmohammadi
%    Omidreza Amrollahinasab
%
% History :
% \change{1.0}{09-Nov-2021}{Original}
%
% --------------------------------------------------
% (c) 2021, Siroos Azizmohammadi,
% Omidreza Amrollahinasab
% Chair of Reservoir Engineering, University of Leoben, Austria
% email: info@dpe.ac.at
% url: dpe.ac.at
% --------------------------------------------------
%
%%
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
        if model.simulation.load_from_sat_prof
            try
                rampup = diff(model.experiment.observation.satProfile.table(2:end,1));
            catch
                fprinf('No saturation profile data provided!')
            end
        else
            rampup = rampupTimesteps(periodInterval, periodInterval/n, rampupsteps);
        end
        schedule = simpleSchedule(rampup, 'bc', model.bc);
        params.scheduleSteps = rampup;
        params.cumScheduleSteps = [params.cumScheduleSteps;...
        cumsum(params.scheduleSteps) + params.cumScheduleSteps(end)];
        schedule.counter = i;
        % simulate the schedule
        %---------------------------------------
%         tic
        verbose = false;
%         [model.state, model.schedulereport] = Simulate(...
%                       model,twoPhaseOilWaterModel,schedule, verbose);
        [~, model.state, model.schedulereport] = evalc...
            ('Simulate(model,twoPhaseOilWaterModel,schedule, verbose)');
%         toc
        %---------------------------------------
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
        if visualization
            params.plotNo = params.plotNo + 1; 
            Visualize(model); 
            frames(params.plotNo) = getframe(gcf); 
        end
        dynamic.params.plotNo = params.plotNo;
    end  
    if visualization
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
%---------------------------------------
% plotting diagnostics toolbar
%     figure
%     plotToolbar(model.grid.G, model.dynamic.states);
end