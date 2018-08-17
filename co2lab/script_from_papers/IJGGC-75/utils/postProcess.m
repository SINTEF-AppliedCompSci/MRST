function postProcess( base, perturbed, varargin )
% Post-processing for parameter uncertainty study


% Calculation of Sørensen-Dice coefficient requires two downloaded
% Mathwork's file exchanges. If not downloaded and available on path, plots
% will be skipped.
polygon_funcs_available = true;
persistent user_warned_once; % if empty, warning can appear, otherwise muted
if ~isdir('Matlab_Polygons_intersection') && ...
        ~isdir('PolygonClipper')
    polygon_funcs_available = false;
    if isempty(user_warned_once)
        user_warned_once = 1;
        fprintf('\n Post-processing of results will skip calculation of\n')
        fprintf(' Sørensen-Dice coefficient due to missing polygon functions.\n\n')
        fprintf(' To avoid this message, ensure you have downloaded the\n')
        fprintf(' following files and added them to your path:\n\n')
        fprintf(' www.mathworks.com/matlabcentral/fileexchange/8818-polygon-clipper\n')
        fprintf(' www.mathworks.com/matlabcentral/fileexchange/18173-polygon_intersection?s_tid=FX_rc2_behav\n\n')
        fprintf(' This message will not be displayed again.\n')
        fprintf(' Continuing after pause short...\n')
        pause(5)
    end
end


% Renaming of INPUTS:
Gt_base         = base.Gt_base;
schedule        = base.schedule;
reports_base    = base.reports_base;
Ma_base         = base.Ma_base;


Gt_all      = perturbed.Gt_all;
rock_all    = perturbed.rock_all;
Ma_all      = perturbed.Ma_all;
reports_all = perturbed.reports_all;

press_deviation_all = [];
if isfield(perturbed,'press_deviation_all')
    press_deviation_all = perturbed.press_deviation_all;
end

opt.surface_pressure    = 1*atm;
opt.NthStep             = 10;
opt.trap_method         = true;
opt.numBins             = [];
opt.year                = []; % can specify year 0 if desired
opt = merge_options(opt, varargin{:});


if ~isempty(opt.year)
    % tstep corresponding to year (tstep = 1 corresponds to year 0, tstep =
    % 2 is first simulated year, etc.)
    years = [0; convertTo(cumsum(schedule.step.val),year)];
    for i=1:numel(opt.year)
        opt.tstep(i) = find(years == opt.year(i));
    end
else
    % the last time step (accounts for initial state/step)
    opt.tstep = numel([reports_all{1}.t]); 
end



% Perturbed models (many realizations):
r = max([ numel(Gt_all), numel(rock_all), numel(press_deviation_all) ]);



if ~isempty(Ma_all)
    
    % ----------------------------------------------------------
    %   forecast curves
    % ----------------------------------------------------------
    fprintf(' Plotting forecast curves...')
    figure, hold on,
    % perturbed
    for i=1:numel(Ma_all)
        plot([0; cumsum(convertTo(schedule.step.val,year))], [Ma_all{i}*1e3], ...
            ':b','LineWidth',1)
    end
    % base
    plot([0; cumsum(convertTo(schedule.step.val,year))], [Ma_base*1e3], ...
        ':b','LineWidth',3)
    ylabel('Amount forecast to stay (Mt)')
    xlabel('Years')
    fprintf(' done.\n')

    
    % ----------------------------------------------------------
    %   plot base trapping inventory with errorbars corresponding to lower
    %       and upper percentiles. Add base forecast curve with errorbars
    % ----------------------------------------------------------
    fprintf(' Plotting trapping inventory...')
    assert( r >= 10, 'Expecting that r >= 10')
    tstep4errorbars = opt.tstep-1; % doesn't account for initial state/step
    tstep4errorbars(tstep4errorbars == 0) = []; % remove possible zeros
    plotTrappingDistribution_withErrorBars(reports_base, Ma_base, ...
        reports_all, Ma_all, schedule, 'lowPercentile',10, ...
        'uppPercentile',90, 'tsteps4errorbars', tstep4errorbars);
    fprintf(' done.\n')


    % ----------------------------------------------------------
    %   plume outlines (perturbed cases and base case)
    % ----------------------------------------------------------
    tstep = opt.tstep;
    if (sum(reports_base(1).sol.s(:,2)) == 0)
        % initial state didn't contain any CO2 so remove any possible
        % tstep=1 from tstep array (i.e., we don't want to make a figure
        % that doesn't contain any plume outlines)
       tstep(tstep == 1) = []; 
    end

    h = figure; hold on;
    % add plume outline in base case (to all subplots)
    plotPlumeOutline(h, Gt_base, reports_base, [tstep], ...
                         'plumeLineWidth', 3, ...
                         'plumeLineStyle', ':', ...
                         'plumeColor', 'r');
    % get base contour of each subplot, one at a time
    hfig = gcf;
    k = numel(tstep):-1:1; % for correct ordering of subplots
    cdata_base = cell(numel(tstep),1);
    for j=1:numel(tstep)
        set(hfig, 'CurrentAxes', hfig.Children(j))
        contours = findobj(hfig.Children(j).Children, 'type','contour');
        % get the latest contour drawn
        cdata_base{k(j),1} = contours(1).ContourMatrix;
        % also add title
        title(['Year ',num2str(opt.year(k(j)))])
    end
    % get contour of set of realizations
    fprintf(' Plotting plume of specified years in realization:\n')
    cdata_realiz = cell(numel(tstep),numel(reports_all));
    for i=1:1:numel(reports_all) % looping thru all realizations can be expensive
        fprintf('  %d of %d\n', i, numel(reports_all))
        plotPlumeOutline(h, Gt_base, reports_all{i}, [tstep], ...
                             'plumeLineWidth', 1, ...
                             'plumeColor', 'b');
        hfig = gcf;
        for j=1:numel(tstep) % loop thru subplots
            set(hfig, 'CurrentAxes', hfig.Children(j))
            contours = findobj(hfig.Children(j).Children, 'type','contour');
            % get latest contour drawn
            cdata_realiz{k(j),i} = contours(1).ContourMatrix;
        end
    end
    % remove empty cells (skipped realizations) and reshape according to
    % number subplots by number contours per subplot
    cdata_realiz = reshape( cdata_realiz(~cellfun('isempty',cdata_realiz)), ...
        numel(tstep), numel(1:1:numel(reports_all)));
    fprintf(' Plume contour data obtained.\n')

    % -------------------------------------------------------------
    %   get SDC from plume outlines (needs to be done before adding grid
    %   outline and base plume outline again)
    % -------------------------------------------------------------
    % loop thru and compute Sørensen-Dice coefficient (SDC) of each
    % obtained contour wrt base contour:
    sdc_realiz = [];
    if polygon_funcs_available
        fprintf(' Calculating Sørensen-Dice coefficient...')
        sdc_realiz = zeros(size(cdata_realiz,1),size(cdata_realiz,2));
        for i=1:size(cdata_realiz,1)
            for j=1:size(cdata_realiz,2)
                [sdc_realiz(i,j), ~] = computeSorensenDiceCoefficient(...
                    cdata_base{i,1}, cdata_realiz{i,j});
            end
        end
        fprintf(' done.\n')
    end
    
    % add base grid outline (to all subplots)
    fprintf(' Adding grid outline...')
    hfig = gcf; hold on;
    for i=1:numel(hfig.Children)
        set(hfig, 'CurrentAxes', hfig.Children(i))
        plotFaces(Gt_base, boundaryFaces(Gt_base))
        axis equal tight
    end
    fprintf(' done.\n')


    % add another base plume outline (on top of realizations)
    fprintf(' Adding base plume outline...')
    plotPlumeOutline(h, Gt_base, reports_base, [tstep], ...
                         'plumeLineWidth', 3, ...
                         'plumeLineStyle', ':', ...
                         'plumeColor', 'r');
    fprintf(' done.\n')

    % --------------------------------------------------------------------
    % Sørensen-Dice coefficient's evolution
    % --------------------------------------------------------------------
    % Plot with average and error bars:
    if ~isempty(sdc_realiz)
        fprintf(' Plotting Sorensen-Dice coefficient...')
        figure, hold on;
        r = size(sdc_realiz,2);
        P10inx = round(r*0.1); if P10inx == 0, P10inx = 1; end;
        P50inx = round(r*0.5);
        P90inx = round(r*0.9);
        [P10, P50, P90, sdc_av, sigma] = deal( zeros(1,numel(tstep)) );
        for i=1:numel(tstep)
            sdc_curr_sorted = sort(sdc_realiz(i,:)); % in ascending order
            P10(i) = sdc_curr_sorted(P10inx);
            P50(i) = sdc_curr_sorted(P50inx);
            P90(i) = sdc_curr_sorted(P90inx);
            sdc_av(i) = mean(sdc_curr_sorted);
            sigma(i) = std(sdc_realiz(i,:));
        end
        if numel(tstep)>1
            myShadedErrorBar(years(tstep), sdc_av, sigma);
        else
            plot(1);
        end
        ylabel('Sørensen-Dice coefficient')
        xlabel('Year')
        fprintf(' done.\n')
    end
    
end


end

function H = myShadedErrorBar(x, y, err)
% modified version of shadedErrorBar.m by Rob Campbell (2009)
% https://se.mathworks.com/matlabcentral/fileexchange/26311-raacampbell-shadederrorbar

    
    H.mainLine = plot(x,y,'-k');
    
    % for patch color
    col = get(H.mainLine,'color');
    edgeColor = col + (1-col)*0.55;
    patchSaturation = 0.15;
    patchColor = col+(1-col)*(1-patchSaturation);
    set(gcf,'renderer','painters')
    
    % plot the patch
    errBar = repmat(err(:)',2,1);
    y = y(:)';
    x = x(:)';
    uE = y + errBar(1,:);
    lE = y - errBar(2,:);
    yP = [lE, fliplr(uE)];
    xP = [x, fliplr(x)];
    H.patch = patch(xP, yP, 1, ...
              'facecolor',patchColor,...
              'edgecolor','none');

    % plot edges around patch
    H.edge(1) = plot(x,lE,'-','color',edgeColor);
    H.edge(2) = plot(x,uE,'-','color',edgeColor);
    
    % re-plot main line since patch covers it up
    H.mainLine = plot(x,y,'-k');

end