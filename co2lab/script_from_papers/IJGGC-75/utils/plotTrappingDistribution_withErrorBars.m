function plotTrappingDistribution_withErrorBars(reports_base, Ma_base, reports_all, Ma_all, schedule, varargin)


    % get reports.masses
    opt.lowPercentile = 10;
    opt.uppPercentile = 90;
    opt.tsteps4errorbars = [];

    opt = merge_options(opt, varargin{:});

    assert(numel(schedule.step.val) >= max(opt.tsteps4errorbars))
    % in this function, tstep = 1 corresponds to first simulated step, not
    % initial state.


    % --------------------------------------
    % Plot trapping inventory of base case
    % --------------------------------------
    h = figure; plot(1); ax = get(h,'currentaxes');
    plotTrappingDistribution(ax, reports_base, 'legend_location', 'southeast');
    fsize = 20;
    set(get(gca, 'xlabel'), 'fontsize', fsize)
    set(get(gca, 'ylabel'), 'fontsize', fsize)
    set(gca,'fontsize', fsize);
    hold on;

    % plot forecast curve
    plot([0; cumsum(convertTo(schedule.step.val,year))], [Ma_base*1e3], ...
        ':b','LineWidth',3, 'DisplayName','Trapping forecast')
    legend('off'); legend('show'); % to refresh legend


    % --------------------------------------
    % Plot error bars according to lower and upper percentile values of
    % trapping wedges based on all realizations, taken at specific time steps
    % --------------------------------------

    r = numel(reports_all); % number of realization cases

    % get P10 and P90 (or P25 and P75) mass values at specific time steps, for
    % each type of trapping mechanism
    inx_low = round(r*opt.lowPercentile/100);
    inx_upp = round(r*opt.uppPercentile/100);


    wedge_vals = []; mass_vals = [];
    Plow_wedge = []; Pupp_wedge = [];
    ttime = [];
    timeStep = [0; cumsum(convertTo(schedule.step.val,year))];
    masses_base = reshape([reports_base.masses]', 8, [])';
    for step = opt.tsteps4errorbars

        [diss, resInStr, res, resInPlume, str, subStr, free, leak] = deal([]);
        masses_all_cum = [];

        for n = 1:numel(reports_all) % loop thru all realizations

            % masses of a realization
            masses_curr = reshape([reports_all{n}.masses]', 8, [])';

            diss        = [diss; masses_curr(step,1)];      % dissolved
            resInStr    = [resInStr; masses_curr(step,2)];  % structural residual
            res         = [res; masses_curr(step,3)];       % residual
            resInPlume  = [resInPlume; masses_curr(step,4)];% residual in plume
            str         = [str; masses_curr(step,5)];       % structural plume
            subStr      = [subStr; masses_curr(step,6)];    % structural subscale
            free        = [free; masses_curr(step,7)];      % free plume
            leak        = [leak; masses_curr(step,8)];      % exited

            masses_all_cum = [masses_all_cum; cumsum(masses_curr(step,:))];
            % all realization's wedge values at given step
        end

        % 
        sorted_diss = sort(diss);
        sorted_resInStr = sort(resInStr);
        sorted_res = sort(res);
        sorted_resInPlume = sort(resInPlume);
        sorted_str = sort(str);
        sorted_subStr = sort(subStr);
        sorted_free = sort(free);
        sorted_leak = sort(leak);


        % wedge and mass value of base case
        wedge_vals = [wedge_vals; cumsum(masses_base(step,:))]; % accumulated
        mass_vals = [mass_vals; masses_base(step,:)]; % not accumulated

        ttime = [ttime; timeStep(step)];


        % lower and upper percentile values for wedges for all realizations (at the curr step):
        sorted_wedges = sort(masses_all_cum,1); % sorts values within columns
        inx = inx_low;
        Plow_wedge = [Plow_wedge; sorted_wedges(inx,1), ...
                        sorted_wedges(inx,2), sorted_wedges(inx,3), sorted_wedges(inx,4), ...
                        sorted_wedges(inx,5), sorted_wedges(inx,6), sorted_wedges(inx,7), ...
                        sorted_wedges(inx,8)];
        inx = inx_upp;
        Pupp_wedge = [Pupp_wedge; sorted_wedges(inx,1), ...
                        sorted_wedges(inx,2), sorted_wedges(inx,3), sorted_wedges(inx,4), ...
                        sorted_wedges(inx,5), sorted_wedges(inx,6), sorted_wedges(inx,7), ...
                        sorted_wedges(inx,8)];

    end
    % plot wedge values with their error bars
    %  NB: 5 and 6th columns get switched in plotTrappingDistribution, which
    %  isn't a problem if subscale trapping is off
    %figure, 
    hold on
    inx = [2 3 4 5 7]; % @@ if dissolved wedge exists, include inx=1
    for i=1:numel(inx)
        inxw = inx(i);
        errorbar(ttime, wedge_vals(:,inxw)./1e9, ...
            (wedge_vals(:,inxw) - Plow_wedge(:,inxw))./1e9, ...
            (Pupp_wedge(:,inxw) - wedge_vals(:,inxw))./1e9, ...
            '.', 'Color','k')
    end


    % --------------------------------------
    % Add error bars to forecast curve according to lower and upper percentiles of forecasts
    % --------------------------------------

    % lower and upper percentiles of forecast for all realizations, taken at
    % specific time steps
    Ma_all_mat = cell2mat(Ma_all);
    FP_low = [];
    FP_upp = [];
    F_base = [];
    for step = opt.tsteps4errorbars %1:NthStep:numel(schedule.step.val)
        %Ma_all_mat(step,:) % forecast values at time step "step" of all realizations
        sorted_Ma = sort(Ma_all_mat(step,:)); % in Gt
        FP_low = [FP_low; sorted_Ma(inx_low)*1000]; % in Mt
        FP_upp = [FP_upp; sorted_Ma(inx_upp)*1000];
        F_base = [F_base; Ma_base(step)*1000];
    end

    % plotting
    hold on;
    errorbar(ttime, F_base, F_base - FP_low, FP_upp - F_base, '.','Color','k');


end





