function main()
%close all;
gravity on;
CO2   = CO2props('rho_big_trunc', ''); 
%CO2 = CO2props('rho_huge', ''); % use this normally
%CO2   = CO2props('rho_small', ''); 

%% Define all necessary parameters for a test case (including grid, fluids,
%  initial conditions, etc.)
incomp_case = setupScenarioInject(CO2, 'compressible' , 'incompressible');
hcomp_case  = setupScenarioInject(CO2, 'compressible' , 'horizontal'    );
comp_case   = setupScenarioInject(CO2, 'compressible' , 'full'          );
% incomp_case = setupScenarioSlope(CO2, 'compressible' , 'incompressible');
% hcomp_case  = setupScenarioSlope(CO2, 'compressible' , 'horizontal'    );
% comp_case   = setupScenarioSlope(CO2, 'compressible' , 'full'          );
% incomp_case = setupScenarioDome(CO2, 'compressible' , 'incompressible');
% hcomp_case  = setupScenarioDome(CO2, 'compressible' , 'incompressible');
% comp_case   = setupScenarioDome(CO2, 'compressible' , 'incompressible');
%% Initialize system structure
icomp_sys = prepareSimulationSystem(incomp_case);
hcomp_sys = prepareSimulationSystem(hcomp_case);    
comp_sys  = prepareSimulationSystem(comp_case);

%% Solving and plotting system

% initial values
incomp_x.pressure    = incomp_case.initial.p0;
incomp_x.h           = incomp_case.initial.h0;

hcomp_x.pressure    = hcomp_case.initial.p0;
hcomp_x.h           = hcomp_case.initial.h0;

comp_x.pressure    = comp_case.initial.p0;
comp_x.h           = comp_case.initial.h0;

% some plot-related variables
xc = incomp_case.Gt.cells.centroids(:,1); % should be same for 'comp_case'
zt = incomp_case.Gt.cells.z;

%keep_timesteps = [5 55 200];
%keep_timesteps = [5 20 40];
%keep_timesteps = [5, 25, 60];
%keep_timesteps = [5, 25, 50];
%keep_timesteps =  [5, 10, 25];
%keep_timesteps = [5, 25, 70];
%keep_timesteps = [5, 115, 220];
keep_timesteps = [21, 60];
saved_timesteps = [];

%% Main time loop
tic;
for t = 1:incomp_case.numSteps
    fprintf('timestep %d\n', t);

    % make changes to well rates/pressures based on the schedule
    incomp_case.Wt = updateWellState(t, incomp_case.schedule, incomp_case.Wt);
    hcomp_case.Wt  = updateWellState(t, hcomp_case.schedule,  hcomp_case.Wt);
    comp_case.Wt   = updateWellState(t, comp_case.schedule,   comp_case.Wt);
    
    [incomp_x , its] = solvefiADI(incomp_x, incomp_case.dt, incomp_case.Wt, incomp_case.Gt, icomp_sys);
    [hcomp_x  , its] = solvefiADI(hcomp_x , hcomp_case.dt , hcomp_case.Wt , hcomp_case.Gt , hcomp_sys);
    [comp_x   , its] = solvefiADI(comp_x  , comp_case.dt  , comp_case.Wt  , comp_case.Gt  ,  comp_sys);
    
    figure(1); locusPlot([comp_x], {'r'}, ...
                         comp_case.gravity.theta, comp_case.temp_grad, ...
                         true, false);    
    figure(2); basicMultiplot(xc, zt, [incomp_x, hcomp_x, comp_x], {'b', 'g', 'r'}, true, 10);%90);    
        
    % textual reporting
    basicReporting('Incompressible case:\n', incomp_case, incomp_x, t == 1, 'incomp');
    basicReporting('Horizontal comp. case:\n', hcomp_case, hcomp_x, t == 1, 'hcomp');
    basicReporting('Compressible case:\n',     comp_case,   comp_x, t == 1, 'comp');
    %pause;
    
    % keeping selected timesteps
    if ismember(t, keep_timesteps)
        saved_timesteps = [saved_timesteps; [incomp_x, hcomp_x, comp_x]];
    end;
end
toc

save('cached', 'saved_timesteps', 'xc', 'zt');

% Plotting graphs for selected timesteps
close all;
num_steps = numel(keep_timesteps);

figure(num_steps + 1); plotsummary(xc, zt, saved_timesteps, {'b', 'g', 'r'}, ...
                                   {'Year 21', 'Year 60'});
% figure(num_steps + 1); plotsummary(xc, zt, saved_timesteps, {'b', 'g', 'r'}, ...
%                                    {'Year 5', 'Year 55', 'Year 140'});

%plotting loci 
%h = PVTplot([5e6, 11e6], [300, 320], CO2);
%h = PVTplot([3e6, 12e6], [290, 327], CO2); % slope
%h = PVTplot([7e6, 10e6], [303, 315], CO2); %dome
h = PVTplot([6e6, 11e6], [302, 316], CO2); %inject
set(gcf, 'position', [0 0 695 499]);
hold on;
 styles = {'k', 'm', 'c'};
for i = 1:num_steps
    % plotting locus
    itemp  = saved_timesteps(i,3).info.intTemp;
    ipress = saved_timesteps(i,3).info.intPress;
    plot(itemp, ipress, styles{i});
end
if false
    h = figure; % new figure
    hold on;
    EOS.compressible = 'full';
    EOS.rho = CO2.rho;
    EOS.beta = CO2.beta;
    EOS.gamma = CO2.gamma;
    EOS.beta2 = @(p, t) CO2.rhoDPP(p,t)./CO2.rho(p,t);
    EOS.gamma2 = @(p, t) CO2.rhoDTT(p,t)./CO2.rho(p,t);
    EOS.chi = @(p, t) CO2.rhoDPT(p,t)./CO2.rho(p,t);
    max_y = 0;
    for i = 1:num_steps;
        % plotting density variation
        int_temp  = saved_timesteps(i,3).info.intTemp;
        int_press = saved_timesteps(i,3).info.intPress;
        tgrad     = comp_case.temp_grad / 1000;
        gcost     = norm(gravity) * cos(comp_case.gravity.theta);
        [~,~,~,~,eta] = etaIntegrals(EOS, int_press, int_temp, tgrad, gcost);
        yvals = eta(-saved_timesteps(i,3).h);
        if max(yvals) > max_y
            max_y = max(yvals);
        end
        
        plot(xc/1e3, yvals , styles{i});
    end
    xlabel('km', 'FontSize', 20);
    ylabel('\rho_T/\rho_M', 'FontSize', 20);
    axis([min(xc)/1e3, max(xc)/1e3, 1, max_y + (max_y-1)*0.15]);
    set(gca, 'FontSize', 16);
end

CO2.dispose(); 
keyboard;
end
%%                                                                             


function geometric1DMultiplot(Gt, states, titles, first_time)

    persistent yrange;
    persistent simnum;
    if first_time
        simnum = numel(states);
        yrange = cell(simnum, 1);
    end
    
    assert(simnum == numel(states));
    assert(Gt.cartDims(1) == prod(Gt.cartDims)); % make sure length is along x direction
    
    xc = Gt.cells.centroids(:, 1)/1e3; % abscissa
    zt = Gt.cells.z;
    zb = zt + Gt.cells.H;
    for i = 1:simnum
        %% Drawing pressure curves
        subplot(simnum, 2, 2*i-1);
        yrange{i} = polygraph([states(i).pressure, ...
                               states(i).info.intPress, ...
                               states(i).info.botPress]/1e6, ...
                              {'r', 'g', 'b'}, ...
                              {'km', 'MPa'}, ...
                              'Top, interface and bottom pressure.', ...
                              xc, ...
                              yrange{i});
        
        %% Drawing grid geometry and height profile
        subplot(simnum, 2, 2*i); cla; hold on;
        hdiff = max(zt) - min(zt) + 5;
        % Drawing caprock
        patch(xc([1 1:end end]), [min(zt)-5; zt; min(zt)-5],.7*[1 1 1]);

        % % Drawing aquifer bottom
        patch(xc([1 1:end end]), [max(zb)+5; zb; max(zb)+5], .7*[1 1 1]);

        % Drawing plume
        z_int = zt + states(i).h;        
        patch(xc([1:end end:-1:1]), [z_int; zt(end:-1:1)], 'r')

        set(gca, 'YDir', 'reverse'), axis tight;
        xlabel('km'); 
        ylabel('depth (m)');
        title(titles{i});
    end
    drawnow;
end
%%                                                                              

function basicMultiplot(abscissa, z_top, states, colors, first_time, H)
    persistent yrange_1;
    persistent yrange_2;
    persistent yrange_3;
    
    if first_time
        % (re-)initialize all static variables
        yrange_1 = [];
        yrange_2 = [];
        yrange_3 = [];
    end
    
    stateinfo = [states.info];
    
    %% plotting top pressure
    subplot(3,4,1);
    yrange_1 = polygraph([stateinfo.topPress]/1e6, ...
                            colors, ...
                            {'km', 'MPa'}, ...
                            'Top pressure', ...
                            abscissa /1e3, ...
                            yrange_1);
    
    %% plotting interface pressure
    subplot(3,4,2);
    yrange_2 = polygraph([stateinfo.intPress]/1e6, ...
                            colors, ...
                            {'km', 'MPa'}, ...
                            'Interface pressure', ...
                            abscissa /1e3, ...
                            yrange_2);
    
    %% plotting CO2 flux
    subplot(3,4,3); yrange_3 = polygraph([stateinfo.fluxCO2], colors, {'km','??'}, ...
                                         'Mass flux', abscissa(2:end)/1e3, yrange_3);
    %% plotting height profile
    subplot(3,4,4);
    heights = bsxfun(@plus, [states.h], z_top);
    max_z = max(z_top);
    min_z = min(z_top);
    polygraph([heights, z_top], {colors{:}, 'k'}, ...
              {'km','m'}, 'Height profile', abscissa/1e3, [min_z, max_z + H]);
    set(gca, 'YDir', 'reverse');
    
    %% plotting interface temperature
    subplot(3,4,5); polygraph([stateinfo.intTemp], colors, {'km', 'K'}, 'Iface temp.', abscissa/1e3, []);
    
    %% plotting interface density
    subplot(3,4,6); polygraph([stateinfo.intRho], colors, {'km', 'kg/m3'}, 'Iface density', abscissa/1e3, []);
    
    %% plotting top density
    subplot(3,4,7); polygraph([stateinfo.topRho], colors, {'km', 'kg/m3'}, 'Top density', abscissa/1e3, []);
    
    %% plotting plume vertical rho difference
    subplot(3,4,8);
    for i = 1:size(stateinfo,2)
        rhoI = stateinfo(i).intRho;
        rhoT = stateinfo(i).topRho;
        rhoDiff(:,i) = (rhoT-rhoI)./rhoI * 100;
    end
    polygraph(rhoDiff, colors, {'km', '%'}, 'Density diff', abscissa /1e3, [-3, 3]);
    
    %% Plotting differences in height profile
    subplot(3,4,9);
    tmp = repmat(states(3).h, 1, 2);
    hdiffrel = [states(1).h, states(2).h] - tmp;
    
    polygraph(hdiffrel, colors, {'km','m'}, 'Height diff', abscissa/1e3, []);
    
    %% Plotting difference in mass flux
    [flux1, flux2, flux3] = stateinfo.fluxCO2;

    fdiff1 = flux1-flux3;
    fdiff2 = flux2-flux3;

    subplot(3, 4, 10); 
    polygraph(fdiff1, colors, {'km', '??'}, 'Flux diff', abscissa(2:end)/1e3, []);
    subplot(3, 4, 11);
    polygraph(fdiff2, colors(2), {'km', '??'}, 'Flux diff', abscissa(2:end)/1e3, []);
    drawnow;
end
%%                                                                              

function basicReporting(intro, simcase, state, first_time, cname)

    persistent injected;
    persistent escaped;
    gct = norm(gravity) * cos(simcase.gravity.theta);     
    Gct = simcase.temp_grad * cos(simcase.gravity.theta)/1000;
    
    if first_time
        % (re-)initialize all static variables
        p0 = simcase.initial.p0;
        h0 = simcase.initial.h0;
        t0 = simcase.top_temp + (Gct .* h0);
        Ieta0 = etaIntegrals(simcase.fluid.CO2, p0, t0, Gct, gct);
        rho_int_0  = simcase.fluid.CO2.rho(p0, t0);
        co2_cols_0 = rho_int_0 .* simcase.initial.h0 .* Ieta0(-h0);
        injected.(cname) = sum(co2_cols_0 .* simcase.Gt.cells.volumes .* simcase.rock.poro);
        escaped.(cname)  = 0;
    end

    % Adding injected CO2 
    rate_ix = find(strcmpi(simcase.Wt.type, 'rate'));
    injected.(cname) = injected.(cname) + simcase.dt * sum(simcase.Wt.val(rate_ix));

    % Registering escaped CO2
    escaped.(cname) = escaped.(cname) + simcase.dt * sum(state.info.outflowCO2);
    
    % Computing contained CO2
    
    T_iface     = simcase.top_temp + (Gct .* state.h);
    T_bot       = simcase.top_temp + (Gct .* simcase.Gt.cells.H);
    Ieta        = etaIntegrals(simcase.fluid.CO2, state.info.intPress, T_iface, Gct, gct);
    rhoInt      = simcase.fluid.CO2.rho(state.info.intPress, T_iface);
    co2_columns = rhoInt .* state.h .* Ieta(-state.h);
    contained   = sum(co2_columns .* simcase.Gt.cells.volumes .* simcase.rock.poro);
    pdiff_TI    = state.info.intPress - state.info.topPress;
    
    % Reporting
    fprintf(intro); fprintf('\n');
    fprintf('Injected so far: %e\n', injected.(cname));
    fprintf('Escaped so far:  %e\n', escaped.(cname));
    fprintf('Left in grid:    %e\n', contained);
    fprintf('Balance:         %e\n', injected.(cname) - escaped.(cname) - contained);
    fprintf('|\n');
    fprintf('Max top press:             %e\n', max(state.info.topPress(state.h > 0)));
    fprintf('Max int press:             %e\n', max(state.info.intPress(state.h > 0)));
    fprintf('Max top-int pressure diff: %e\n', max(pdiff_TI(state.h > 0)))
    fprintf('Max bot press:             %e\n', max(state.info.botPress(state.h > 0)));
    fprintf('|\n');
    if (strcmpi(simcase.fluid.CO2.compressible, 'full'))
        rhoTop = simcase.fluid.CO2.rho(state.info.topPress, simcase.top_temp);
    else        
        rhoTop = rhoInt;
    end
    fprintf('Max rho on top: %e\n', max(rhoTop(state.h > 0)));
    fprintf('Max rho on int: %e\n', max(rhoInt(state.h > 0)));
    fprintf('Av. plume rho:  %e\n', sum(co2_columns) / sum(state.h));
    % bRho = simcase.fluid.CO2.rho(state.info.botPress, T_bot);
    % fprintf('Mean rho on bot: %e\n', mean(bRho(state.h > 0)));
    fprintf('$\n');
    fprintf('-------------------------\n\n');
    
end
%%                                                                              


% function yrange = polygraph(graphs, colors, labels, plot_title, xvals, yrange)

%     hold on; cla;
%     for g_ix = 1:size(graphs, 2) % one graph per column
        
%         plot(xvals, graphs(:, g_ix), colors{g_ix});
        
%     end
%     xlabel(labels{1});
%     ylabel(labels{2});

%     if isempty(yrange) || max(max(graphs)) > yrange(2) || min(min(graphs)) < yrange(1)
%         % Y-range not yet set/imposed.  Computing a reasonable range
%         PADDING = 0.1;
        
%         yrange = [min(min(graphs)), max(max(graphs))];
%         span = yrange(2) - yrange(1);
%         yrange(1) = yrange(1) - PADDING * span;
%         yrange(2) = yrange(2) + PADDING * span;
    
%     end        
   
%     title(plot_title);
    
%     % small 'hack' to avoid collapse if y is constant (range of zero...)
%     if yrange(1) < yrange(2)
%         margin = 0;
%     else
%         margin = yrange(1) * 0.05;
%     end
%     axis([min(xvals), max(xvals), yrange(1)-margin, yrange(2)+margin]);
    
% end
%%                                                                              

% update the well state based on timestep and schedule
function Wt = updateWellState(tstep, schedule, Wt)
    
    now_ix        = find(schedule(:,1) == tstep);
    now_changes   = schedule(now_ix, 2:3);
    wells         = now_changes(:,1);
    vals          = now_changes(:,2);

    Wt.val(wells) = vals;
    
end
%%                                                                              

function tcase = setupScenario1(CO2, varargin)
%% Injection into horizontal aquifer
    
    % Setting basic parameters
    temp_grad = 55;        
    
        
    opt = struct('cartDims',      [50 1 1], ...
                 'porosity',      0.1, ...
                 'perm',          40 * milli * darcy, ... 
                 'physDims',      [5000 3000 100], ...
                 'topo',          1000, ...
                 'gravityTheta',  (0*pi/180), ...
                 'wellPos',       [25, 1], ...
                 'h0',            0, ...
                 'ref_depth',     0, ...
                 'ref_press',     1*atm, ...
                 'ref_temp',      273.15 + 4 - 0.1 * temp_grad, ...
                 'ref_cell',      [1, 1], ...
                 'temp_grad',     temp_grad, ... % deg/km 
                 'mu',            [5.36108e-5, 6.5e-4], ...
                 'waterDensity',  1200);
    
    % Adding cell structure parameters
    opt.bcondTypes   = {'Pressure', 'Pressure', 'Flux', 'Flux'};
    opt.wellType     = {'Rate'};
    opt.compressible = []; % no default - provide by 'varargin'
    % opt.schedule     = {[1, 0.4e6 * kilo * kilogram / year; ...
    %                      10, 0]};
    opt.schedule     = {[1, 1e6 * kilo * kilogram / year; ...
                         2, 0]};
    opt.simTime      = 50 * year;
    opt.timesteps    = 50; % 3-month long timesteps

    opt = merge_options(opt, varargin{:});
    cell_opt = fullStructToCell(opt);
    tcase = setupSharpInterfaceCartesianTestCase(CO2, cell_opt{:});
end
%%                                                                              

function tcase = setupScenarioInject(CO2, varargin)
%% Injection into a horizontal, plane aquifer
    resolution = 100;
    % Setting basic parameters
    temp_grad = 40;     
    
    h0 = zeros(resolution, 1);

    opt = struct('cartDims',      [resolution 1 1], ...
                 'porosity',      0.1, ...
                 'perm',          400 * milli * darcy, ...  % 400
                 'physDims',      [40000 3000 150], ...
                 'topo',          750, ...
                 'gravityTheta',  (0*pi/180), ...
                 'wellPos',       [resolution/2, 1], ...
                 'h0',            h0, ...
                 'ref_depth',     0, ...
                 'ref_press',     1*atm, ...
                 'ref_temp',      273.15 + 6, ...% - 0.1 * temp_grad, ...
                 'ref_cell',      [resolution/2, 1], ...
                 'temp_grad',     temp_grad, ... % deg/km 
                 'mu',            [5.36108e-5, 5.4e-5], ... %6.5e-4], ...
                 'waterDensity',  1050);
    
    % Adding cell structure parameters
    opt.bcondTypes   = {'Pressure', 'Pressure', 'Flux', 'Flux'};
    opt.wellType     = {'Rate'};
    opt.compressible = []; % no default - provide by 'varargin'
    %opt.schedule     = {[1, 0]};
    opt.schedule     = {[1, 1e7 * kilo * kilogram / year; 20, 0]};
    opt.simTime      = 60 * year;
    opt.timesteps    = 60;

    opt = merge_options(opt, varargin{:});
    cell_opt = fullStructToCell(opt);
    tcase = setupSharpInterfaceCartesianTestCase(CO2, cell_opt{:});
end

function tcase = setupScenarioInject_orig(CO2, varargin)
%% Injection into a horizontal, plane aquifer
    resolution = 100;
    % Setting basic parameters
    temp_grad = 45;     
    
    h0 = zeros(resolution, 1);
    
        
    opt = struct('cartDims',      [resolution 1 1], ...
                 'porosity',      0.1, ...
                 'perm',          500 * milli * darcy, ... 
                 'physDims',      [40000 3000 150], ...
                 'topo',          800, ...
                 'gravityTheta',  (0*pi/180), ...
                 'wellPos',       [resolution/2, 1], ...
                 'h0',            h0, ...
                 'ref_depth',     0, ...
                 'ref_press',     1*atm, ...
                 'ref_temp',      273.15 + 6, ...% - 0.1 * temp_grad, ...
                 'ref_cell',      [resolution/2, 1], ...
                 'temp_grad',     temp_grad, ... % deg/km 
                 'mu',            [5.36108e-5, 5.4e-5], ... %6.5e-4], ...
                 'waterDensity',  1100);
    
    % Adding cell structure parameters
    opt.bcondTypes   = {'Pressure', 'Pressure', 'Flux', 'Flux'};
    opt.wellType     = {'Rate'};
    opt.compressible = []; % no default - provide by 'varargin'
    %opt.schedule     = {[1, 0]};
    opt.schedule     = {[1, 1e7 * kilo * kilogram / year; 10, 0]};
    opt.simTime      = 80 * year;
    opt.timesteps    = 60;

    opt = merge_options(opt, varargin{:});
    cell_opt = fullStructToCell(opt);
    tcase = setupSharpInterfaceCartesianTestCase(CO2, cell_opt{:});
end


function tcase = setupScenarioInjectKrull(CO2, varargin)
%% Injection into a horizontal, plane aquifer
    resolution = 100;
    % Setting basic parameters
    temp_grad = 0;     
    
    h0 = zeros(resolution, 1);
    
        
    opt = struct('cartDims',      [resolution 1 1], ...
                 'porosity',      0.1, ...
                 'perm',          500 * milli * darcy, ... 
                 'physDims',      [40000 3000 150], ...
                 'topo',          800, ...
                 'gravityTheta',  0, ...
                 'wellPos',       [resolution/2, 1], ...
                 'h0',            h0, ...
                 'ref_depth',     0, ...
                 'ref_press',     1*atm, ...
                 'ref_temp',      273.15 + 35, ...% - 0.1 * temp_grad, ...
                 'ref_cell',      [resolution/2, 1], ...
                 'temp_grad',     temp_grad, ... % deg/km 
                 'mu',            [5.36108e-5, 5.4e-5], ... %6.5e-4], ...
                 'waterDensity',  1200);
    
    % Adding cell structure parameters
    opt.bcondTypes   = {'Pressure', 'Pressure', 'Flux', 'Flux'};
    opt.wellType     = {'Rate'};
    opt.compressible = []; % no default - provide by 'varargin'
    %opt.schedule     = {[1, 0]};
    opt.schedule     = {[1, 1e7 * kilo * kilogram / year; 20, 0]};
    opt.simTime      = 80 * year;
    opt.timesteps    = 60;

    opt = merge_options(opt, varargin{:});
    cell_opt = fullStructToCell(opt);
    tcase = setupSharpInterfaceCartesianTestCase(CO2, cell_opt{:});
end



% function tcase = setupScenarioSlope(CO2, varargin)
% %% Tall plume flowing up a sloping aquifer
%     resolution = 200;%200;%100;
%     % Setting basic parameters
%     temp_grad = 45;%22;     
%     h0 = zeros(resolution, 1);
    
%     h0(20 : 40 ) = - linspace(0, 100, 21) * 4e7;
%     h0(40 : 60) = - linspace(100, 0, 21) * 4e7;

%     %h0(resolution/5 : 3*resolution/10) = 100;
        
%     opt = struct('cartDims',      [resolution 1 1], ...
%                  'porosity',      0.1, ...
%                  'perm',          1000 * milli * darcy, ... 
%                  'physDims',      [80000 3000 150], ... %[40000 3000 150], ...%[20000 3000 150], ...
%                  'topo',          800, ...
%                  'gravityTheta',  ((1/2)*pi/180), ...
%                  'wellPos',       [resolution/2, 1], ...
%                  'h0',            h0, ...
%                  'ref_depth',     0, ...
%                  'ref_press',     1*atm, ...
%                  'ref_temp',      273.15 + 6,...%17, ...% - 0.1 * temp_grad, ...
%                  'ref_cell',      [1, 1], ...
%                  'temp_grad',     temp_grad, ... % deg/km 
%                  'mu',            [5.36108e-5, 6.5e-4], ...
%                  'waterDensity',  1100);
%     % Adding cell structure parameters
%     opt.bcondTypes   = {'Pressure', 'Pressure', 'Flux', 'Flux'};
%     opt.wellType     = {'Rate'};
%     opt.compressible = []; % no default - provide by 'varargin'
%     opt.schedule     = {[1, 0]};
%     %opt.schedule     = {[1, 1e8 * kilo * kilogram / year; 2, 0]};
%     opt.simTime      = 50 * year;
%     opt.timesteps    = 50;

%     opt = merge_options(opt, varargin{:});
%     cell_opt = fullStructToCell(opt);
%     tcase = setupSharpInterfaceCartesianTestCase(CO2, cell_opt{:});
% end



function tcase = setupScenarioSlope(CO2, varargin)
%% Tall plume flowing up a sloping aquifer
    resolution = 200;%200;%100;
    % Setting basic parameters
    temp_grad = 45;
    h0 = zeros(resolution, 1);
    
    h0(20 : 40 ) = - linspace(0, 100, 21) * 4e7;
    h0(40 : 60) = - linspace(100, 0, 21) * 4e7;

    %h0(resolution/5 : 3*resolution/10) = 100;
        
    opt = struct('cartDims',      [resolution 1 1], ...
                 'porosity',      0.1, ...
                 'perm',          1400 * milli * darcy, ... 
                 'physDims',      [80000 3000 150], ... %[40000 3000 150], ...%[20000 3000 150], ...
                 'topo',          1000, ...
                 'gravityTheta',  ((1/2)*pi/180), ...
                 'wellPos',       [resolution/2, 1], ...
                 'h0',            h0, ...
                 'ref_depth',     0, ...
                 'ref_press',     1*atm, ...
                 'ref_temp',      273.15 + 8,...%17, ...% - 0.1 * temp_grad, ...
                 'ref_cell',      [1, 1], ...
                 'temp_grad',     temp_grad, ... % deg/km 
                 'mu',            [5.36108e-5, 6.5e-4], ...
                 'waterDensity',  1100);
    % Adding cell structure parameters
    opt.bcondTypes   = {'Pressure', 'Pressure', 'Flux', 'Flux'};
    opt.wellType     = {'Rate'};
    opt.compressible = []; % no default - provide by 'varargin'
    opt.schedule     = {[1, 0]};
    %opt.schedule     = {[1, 1e8 * kilo * kilogram / year; 2, 0]};
    opt.simTime      = 50 * year;
    opt.timesteps    = 50;

    opt = merge_options(opt, varargin{:});
    cell_opt = fullStructToCell(opt);
    tcase = setupSharpInterfaceCartesianTestCase(CO2, cell_opt{:});
end




% function tcase = setupScenarioSlope(CO2, varargin)
% %% Tall plume flowing up a sloping aquifer
%     resolution = 250;
%     % Setting basic parameters
%     temp_grad = 40;%22;     
%     h0 = zeros(resolution, 1);
    
%     h0(40 : 50 ) = - linspace(0, 130, 11) * 4e7;
%     h0(50 : 60)  = - linspace(130, 0, 11) * 4e7;

%     %h0(resolution/5 : 3*resolution/10) = 100;
        
%     opt = struct('cartDims',      [resolution 1 1], ...
%                  'porosity',      0.1, ...
%                  'perm',          1000 * milli * darcy, ... 
%                  'physDims',      [50000 3000 150], ... %[40000 3000 150], ...%[20000 3000 150], ...
%                  'topo',          800, ...
%                  'gravityTheta',  ((1/3)*pi/180), ...
%                  'wellPos',       [resolution/2, 1], ...
%                  'h0',            h0, ...
%                  'ref_depth',     0, ...
%                  'ref_press',     1*atm, ...
%                  'ref_temp',      273.15 + 6,...%17, ...% - 0.1 * temp_grad, ...
%                  'ref_cell',      [1, 1], ...
%                  'temp_grad',     temp_grad, ... % deg/km 
%                  'mu',            [5.36108e-5, 6.5e-4], ...
%                  'waterDensity',  1100);
%     % Adding cell structure parameters
%     opt.bcondTypes   = {'Pressure', 'Pressure', 'Flux', 'Flux'};
%     opt.wellType     = {'Rate'};
%     opt.compressible = []; % no default - provide by 'varargin'
%     opt.schedule     = {[1, 0]};
%     %opt.schedule     = {[1, 1e8 * kilo * kilogram / year; 2, 0]};
%     opt.simTime      = 70 * year;
%     opt.timesteps    = 70;

%     opt = merge_options(opt, varargin{:});
%     cell_opt = fullStructToCell(opt);
%     tcase = setupSharpInterfaceCartesianTestCase(CO2, cell_opt{:});
% end

% function tcase = setupScenarioSlope(CO2, varargin)
% %% Tall plume flowing up a sloping aquifer
%     resolution = 200;
%     % Setting basic parameters
%     temp_grad = 35;%22;     
%     h0 = zeros(resolution, 1);
    
%     h0(30 : 40 ) = - linspace(0, 100, 11) * 4e7;
%     h0(40 : 50)  = - linspace(100, 0, 11) * 4e7;

%     %h0(resolution/5 : 3*resolution/10) = 100;
        
%     opt = struct('cartDims',      [resolution 1 1], ...
%                  'porosity',      0.1, ...
%                  'perm',          1000 * milli * darcy, ... 
%                  'physDims',      [100000 3000 150], ... %[40000 3000 150], ...%[20000 3000 150], ...
%                  'topo',          700, ...
%                  'gravityTheta',  ((1/3)*pi/180), ...
%                  'wellPos',       [resolution/2, 1], ...
%                  'h0',            h0, ...
%                  'ref_depth',     0, ...
%                  'ref_press',     1*atm, ...
%                  'ref_temp',      273.15 + 12,...%17, ...% - 0.1 * temp_grad, ...
%                  'ref_cell',      [1, 1], ...
%                  'temp_grad',     temp_grad, ... % deg/km 
%                  'mu',            [5.36108e-5, 6.5e-4], ...
%                  'waterDensity',  1200);
%     % Adding cell structure parameters
%     opt.bcondTypes   = {'Pressure', 'Pressure', 'Flux', 'Flux'};
%     opt.wellType     = {'Rate'};
%     opt.compressible = []; % no default - provide by 'varargin'
%     opt.schedule     = {[1, 0]};
%     %opt.schedule     = {[1, 1e8 * kilo * kilogram / year; 2, 0]};
%     opt.simTime      = 25 * year;
%     opt.timesteps    = 25;

%     opt = merge_options(opt, varargin{:});
%     cell_opt = fullStructToCell(opt);
%     tcase = setupSharpInterfaceCartesianTestCase(CO2, cell_opt{:});
% end

%%                                                                              


% function tcase = setupScenarioDome(CO2, varargin)
% %% Injection into dome, close to critical point
%     % Setting basic parameters
%     resolution = 50;
%     top_geom = +60 * sin([1:(resolution+1)]'/(resolution+1) * 3 * pi) + 800;
%     % top_geom(top_geom < 0) = -30;
%     % top_geom(top_geom > 0) = 30;
%     % top_geom(23:28) = -60;
    
%     sample_h_profile = zeros(50, 1); h_profile(25:35) = 30;
    
%     opt = struct('cartDims',      [resolution 1 1], ...
%                  'porosity',      0.1, ...
%                  'perm',          2000 * milli * darcy, ...
%                  'physDims',      [10000 3000 100], ...
%                  'topo',          [top_geom, top_geom], ...
%                  'gravityTheta',  (0*pi/180), ...
%                  'wellPos',       [ceil(resolution/2), 1], ...
%                  'h0',            0, ...%sample_h_profile, ... %0, ...
%                  'ref_depth',     0, ...
%                  'ref_press',     1*atm, ...
%                  'ref_temp' ,     273.15 + 4, ...
%                  'ref_cell' ,     [ceil(resolution/4), 1], ...
%                  'temp_grad',     40, ...
%                  'mu',            [5.36108e-5, 6.5e-4], ...
%                  'waterDensity',  1000);

%     % Adding cell structure parameters
%     opt.bcondTypes   = {'Pressure', 'Pressure', 'Flux', 'Flux'};
%     opt.wellType     = {'Rate'};
%     opt.compressible = []; % no default - provide by 'varargin'
%     opt.schedule     = {[1, 2.5e6 * kilo * kilogram / year; 40 0]};
%     opt.simTime      = 30 * year;
%     opt.timesteps    = 120; %120; % 3-month long timesteps

%     opt = merge_options(opt, varargin{:});
%     cell_opt = fullStructToCell(opt);
%     tcase = setupSharpInterfaceCartesianTestCase(CO2, cell_opt{:});
% end
%%                                                                              


% ============================== THIS ONE WORKS ==============================
% function tcase = setupScenarioDome(CO2, varargin)
% %% Injection into dome, close to critical point
%     % Setting basic parameters
%     resolution = 100;
%     top_geom = +60 * sin([1:(resolution+1)]'/(resolution+1) * 3 * pi) + 800;
%     % top_geom(top_geom < 0) = -30;
%     % top_geom(top_geom > 0) = 30;
%     % top_geom(23:28) = -60;
    
%     opt = struct('cartDims',      [resolution 1 1], ...
%                  'porosity',      0.1, ...
%                  'perm',          1000 * milli * darcy, ...
%                  'physDims',      [10000 3000 100], ...
%                  'topo',          [top_geom, top_geom], ...
%                  'gravityTheta',  (0*pi/180), ...
%                  'wellPos',       [ceil(resolution/2), 1], ...
%                  'h0',            0, ...%sample_h_profile, ... %0, ...
%                  'ref_depth',     0, ...
%                  'ref_press',     1*atm, ...
%                  'ref_temp' ,     273.15 + 4, ...
%                  'ref_cell' ,     [ceil(resolution/4), 1], ...
%                  'temp_grad',     40, ...
%                  'mu',            [5.36108e-5, 6.5e-4], ...
%                  'waterDensity',  1000);

%     % Adding cell structure parameters
%     opt.bcondTypes   = {'Pressure', 'Pressure', 'Flux', 'Flux'};
%     opt.wellType     = {'Rate'};
%     opt.compressible = []; % no default - provide by 'varargin'
%     opt.schedule     = {[1, 1.5e6 * kilo * kilogram / year; 10 0]};
%     opt.simTime      = 40 * year;
%     opt.timesteps    = 40; %120; % 3-month long timesteps

%     opt = merge_options(opt, varargin{:});
%     cell_opt = fullStructToCell(opt);
%     tcase = setupSharpInterfaceCartesianTestCase(CO2, cell_opt{:});
% end
% %%                                                                              

function tcase = setupScenarioDome(CO2, varargin)
%% Injection into dome, close to critical point
    % Setting basic parameters
    resolution = 100;
    top_geom = +80 * cos([1:(resolution+1)]'/(resolution+1) * 2 * pi) + 780;
    % ensuring top-geometry is exactly symmetric around center
    top_geom = (top_geom + fliplr(top_geom))/2;
    
    % top_geom(top_geom < 0) = -30;
    % top_geom(top_geom > 0) = 30;
    % top_geom(23:28) = -60;
    
    opt = struct('cartDims',      [resolution 1 1], ...
                 'porosity',      0.2, ...
                 'perm',          1100 * milli * darcy, ...
                 'physDims',      [10000 3000 100], ...
                 'topo',          [top_geom, top_geom], ...
                 'gravityTheta',  (0*pi/180), ...
                 'wellPos',       [ceil(resolution/2), 1], ...
                 'h0',            0, ...%sample_h_profile, ... %0, ...
                 'ref_depth',     0, ...
                 'ref_press',     1*atm, ...
                 'ref_temp' ,     273.15 + 6, ...
                 'ref_cell' ,     [ceil(resolution/4), 1], ...
                 'temp_grad',     40, ...
                 'mu',            [5.36108e-5, 6.5e-4], ...
                 'waterDensity',  1100);

    % Adding cell structure parameters
    %opt.bcondTypes   = {'Pressure', 'Pressure', 'Flux', 'Flux'};
    opt.bcondTypes   = {'Pressure', 'Pressure', 'Flux', 'Flux'};
    opt.wellType     = {'Rate'};
    opt.compressible = []; % no default - provide by 'varargin'
    opt.schedule     = {[1, 2e6 * kilo * kilogram / year; 50 0]};
    opt.simTime      = 200 * year;
    opt.timesteps    = 200; %120; % 3-month long timesteps

    opt = merge_options(opt, varargin{:});
    cell_opt = fullStructToCell(opt);
    tcase = setupSharpInterfaceCartesianTestCase(CO2, cell_opt{:});
end
%%                                                                              
