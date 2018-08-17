%% Comparison of fault cases simulated in Stø aquifer
% The following script generates results (similar**) to those presented in
% Figures 13, 14, and 15 of the paper (see README.txt).

% The Stø aquifer contains several faults. Fault data was provided to us
% through personal communication with the Norwegian Petroleum Directorate.
% We treat these faults as either conducting, semi-conducting, or sealing.
% An injection scenario with 7 wells operating at the same fixed rate is
% simulated, followed by a migration period of close to 3000 years. Over
% pressure induced by the injection scenario is assessed, as well as plume
% migration.

% ** Note: results presented in the paper were generated using a high
% resolution of the grid and relatively small time step sizes. As such the
% time required to simulate the following 30 cases can take several hours.
% To speed up these computations, the following can be done:
%   * use a coarser grid by setting "coarse_level" to 4 or 5
%   * simulate only 100 years post-injection instead of ~3000 years
% These changes will allow you to obtain results quicker, however they are
% likely to be slightly different from those presented in the paper.

mrstModule add co2lab
moduleCheck('ad-core');

%% Model set-up
% --------------------------------------------
% Set up model with faults
% --------------------------------------------
coarse_level = 4; % @@ set to 1 to match same grid resolution used in paper
                  % @@ otherwise try a value of 4 or 5 to speed up simulations
[Gt, rock2D, faults, faces, faultFaces2D] = ...
    setUpFaultExample('coarse_level',coarse_level);


% --------------------------------------------
% Set up simulation details
% --------------------------------------------
shorter_sim_run = true; % @@ set to false to match same schedule used in paper
                        % @@ otherwise simulate for a shorter migration period
[initState, schedule, fluid, seainfo] = setUpForSimulation(Gt, rock2D, ...
    'shorter_sim_run',shorter_sim_run);

% Set up 7 wells such that they are upstream from a fault.
pt1 = [9.247e5, 7.983e6];
pt2 = [9.307e5, 7.985e6];
pt3 = [9.397e5, 7.988e6];
pt4 = [9.472e5, 7.989e6];
pt5 = [9.577e5, 7.992e6];
pt6 = [9.142e5, 7.937e6];
pt7 = [9.097e5, 7.942e6];
wc = findEnclosingCell( Gt, [pt1; pt2; pt3; pt4; pt5; pt6; pt7] );


%% Simulate all cases
% Three different transmissibility multipliers are considered (to treat
% faults as conducting, semi-sealing, or sealing). Also, 10 different
% injection rates applied to the wells are considered.

% for saving results of 3 transMult cases
myMainFolder = fullfile(mrstOutputDirectory, 'simFaults_Sto'); 

transMults = [1 0.01 0]; % transmissibility multipliers for fault-faces
% 0     - sealing
% 0.01  - semi-sealing / semi-conducting
% 1     - no faults (i.e., conducting)

for k = 1:numel(transMults)
    
    transMult = transMults(k);
    vals = [0.01:0.01:0.1]; % m3/s, injection rate per well
    
    for j=1:numel(vals)
        
        val = vals(j);

        % wells:
        W = [];
        for i=1:numel(wc)
            W = addWell(W, Gt.parent, rock2D, wc(i), ...
                        'name',     ['Winj' num2str(wc(i))],  ...
                        'Type',     'rate', ...
                        'Val',      val, ...
                        'comp_i',   [0,1], ...
                        'Radius',   0.3);
        end
        W = convertwellsVE(W, Gt.parent, Gt, rock2D, 'ip_tpf');
        W_shut = W;
        for i=1:numel(wc)
            W_shut(i).val = sqrt(eps);
        end
        schedule.control(1).W = W;
        schedule.control(2).W = W_shut;


        % Construct model, then modify transmissibilites appropriately.
        model = CO2VEBlackOilTypeModel(Gt, rock2D, fluid);
        T_all                 = getVerticallyAveragedTrans(Gt, rock2D);
        T_all_orig            = T_all;
        T_all( faultFaces2D ) = T_all( faultFaces2D ).*transMult;
        model.operators.T_all = T_all;
        model.operators.T     = T_all( model.operators.internalConn );


        % Simulate.
        [wellSols, states] = simulateScheduleAD(initState, model, schedule);


        % Save results for post-processing
        ta = trapAnalysis(Gt, false);
        reports = makeReports(Gt, {initState, states{:}}, rock2D, fluid, ...
                    schedule, [model.fluid.res_water, model.fluid.res_gas], ...
                    ta, []);

        folder = fullfile(myMainFolder, ['transMult' num2str(transMult)]);
        if ~exist(folder,'dir')
            mkdir(folder);
        end
        filename = fullfile(folder, ['injRate' num2str(val) '.mat']);
        save(filename, 'Gt','reports','schedule','rock2D','faultFaces2D', ...
            '-v7.3')

    end
end


%% Make plots shown in Figure 13 of paper:
% Set-up of faults and wells

figure, hold on

% show grid:
plotGrid(Gt,'edgealpha',0.1)

% show bdrys:
bdrys = schedule.control(1).bc.face;
plotFaces(Gt, bdrys, 'edgec','g','linewidth',3)

% show faults:
% (We also plot fault lines 48 and 39, but note that they were not included
% as faults in simulation because they are part of a parallel fault. See
% note on handling parallel faults in setUpFaultExample.m)
plotFaces(Gt, faultFaces2D, 'edgec','b','linewidth',3)
plotFaces(Gt, faces{48}, 'edgec','b','linewidth',1)
plotFaces(Gt, faces{39}, 'edgec','b','linewidth',1)
axis tight off, daspect([1 1 0.09]), view([-120, 30])
light('Position',[-1 -1 -1], 'Style','infinite'); lighting phong

% show wells:
W = schedule.control(1).W;
[~, htext] = plotWell(Gt.parent, W, 'color','k');
replaceWellLabels = [schedule.control(1).W.cells];
for i=1:numel(htext)
    curr_label = get(htext(i),'String');
    if all(ismember(curr_label(5:end), num2str(replaceWellLabels(i))))
        set(htext(i),'String',['W',num2str(i)], 'HorizontalAlignment','center')
    end
    if i==6
        set(htext(i), 'HorizontalAlignment','left')
    elseif i==7
        set(htext(i), 'HorizontalAlignment','right')
    end
end

% add scale to plot:
axis on
box on
[x0, y0, z0] = deal(950000, 8050000, 3000);
line([x0, x0+10000], [y0, y0], [z0, z0])
line([x0, x0], [y0, y0-10000], [z0, z0])
line([x0, x0], [y0, y0], [z0, z0-1000])
axis off


% Another figure to show faults in Hammerfest Aquifer Basin around Stø
figure, hold on
for f = 1:size(faults,2)
    hl = line(faults{f}(:,1), faults{f}(:,2), 'linewidth',2);
end
hg1 = plotGrid(Gt, 'edgealpha',0.1, 'facec','y');
axis equal tight, box on, grid on
legend([hg1, hl],{'Stø','fault lines'},'Location','SouthEast','FontSize',18)
set(gcf,'Position',[3164 420 614 694])
set(gca,'FontSize',12)
xlim([0.875e6 1e6])
ylim([7.90e6 8.06e6])


%% Make plots shown in Figure 14 & 15 of paper:
% Injection rates vs. max over-pressure reached
% For each of the cases simulated above, we calculate the maximum
% over-pressure simulated and plot the results in a line plot. We also
% illustrate results of the 0.25 Mt/yr injection rate case: (i) map of
% over-pressure at the year when maximum over-pressure was reached, (ii)
% map of CO2 saturations after 3000 years.

hh = 102; % figure handle for line plot

folders = { fullfile(myMainFolder, 'transMult0'); ...
            fullfile(myMainFolder, 'transMult0.01'); ...
            fullfile(myMainFolder, 'transMult1') };

for f=1:numel(folders)
    
    folder = folders{f};
    files = dir(folder);
    
    fprintf(' Plotting any data found in folder "%s" ...', folder)
    
    for j=1:numel(files)

        if (files(j).isdir)
           % do nothing

        else

            % load file
            load(fullfile(folder,files(j).name))
            
            inj_rate = schedule.control(1).W.val; % assume all w/ same rate
            states = {reports.sol}';
            p_init = reports(1).sol.pressure;
            
            % determine how much the overburden pressure was surpassed
            maxOverP = zeros(1,numel(reports));
            loc_ind = zeros(1,numel(reports));
            for i=1:numel(reports)

                % pressure field at i-th time step
                p = reports(i).sol.pressure;

                % over pressure
                overP = p - p_init;

                % highest over pressure reached in domain at this time
                % step, and its location
                [maxOverP(i), loc_ind(i)] = max(overP);

            end

            % highest over pressure reached in domain over all time steps,
            % and its location
            [maxMaxOverP, t_ind] = max(maxOverP);
            cell_ind = loc_ind(t_ind); % not cell number, but cell index

            % over pressure at time step when max over pressure was reached
            maxOverP_field = reports(t_ind).sol.pressure - p_init;
            time_step1 = convertTo(reports(t_ind).t,year);

            
            % compare max overP reached vs rates
            figure(hh), hold on
            [~, name, exten] = fileparts(folder);
            fname = [name,exten];
            if strcmpi(fname,'transMult0')
                h1 = plot(inj_rate * 760/1e9*(365*24*60*60), ...
                    convertTo(maxMaxOverP,barsa), 'xk');
                str = 'Sealing case';
            elseif strcmpi(fname,'transMult0.01')
                h2 = plot(inj_rate * 760/1e9*(365*24*60*60), ...
                    convertTo(maxMaxOverP,barsa), 'or');
                str = 'Semi-sealing case';
            elseif strcmpi(fname,'transMult1')
                h3 = plot(inj_rate * 760/1e9*(365*24*60*60), ...
                    convertTo(maxMaxOverP,barsa), '+b');
                str = 'Conducting case';
            end
            
            % make map plots of over-pressure and CO2 saturations for the
            % case with an injection rate of 0.01 m3/s
            if contains(files(j).name,'injRate0.01')

                % over pressure at time step when max over pressure was
                % reached
                figure, hold on
                plotCellData(Gt, convertTo(maxOverP_field,barsa), ...
                    'edgecolor','none')
                plotCellData(Gt, convertTo(maxOverP_field,barsa), ...
                    cell_ind, 'facecolor','red','edgecolor','none')
                colorbar
                title({str;'Over pressure';['year ',num2str(time_step1)]})
                colormap('hot')
                axis tight equal off

                % saturations at end of simulated period (approx 3000
                % years)
                figure
                plotCellData(Gt, states{end}.s(:,2), ...
                    states{end}.s(:,2) > 0.01, 'edgecolor','none');
                axis equal tight off; colorbar
                plotGrid(Gt, 'facecolor','none', 'edgealpha',0.1)
                if strcmpi(str,'Semi-sealing case') || ...
                        strcmpi(str,'Sealing case')
                    plotFaces(Gt, faultFaces2D, 'edgec','r','linewidth',2)
                end
                title({str;'CO2 saturations'; ...
                    ['year ',num2str(convertTo(reports(end).t,year))]})

            end
            
        end

    end
    
    fprintf(' done\n')
    
end
figure(hh), hold on;
title('x - sealing, o - semi-sealing, + - conducting')
xlabel('Injection rate (Mt/yr)')
ylabel('Max over pressure (bars)')
box on

% add hypothetical over-pressure limit
plot([0 2.5], [70 70], 'k--')