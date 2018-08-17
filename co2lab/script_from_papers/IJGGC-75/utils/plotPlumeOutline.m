function plotPlumeOutline(h, Gt, reports, plot_steps, varargin)
% uses mapPlot to plot a plume outline at specific time steps in separate
% subplots

    opt.plume_threshold     = 0.3;
    opt.plumeLineWidth      = 2;
    opt.plumeLineStyle      = '-';
    opt.plumeColor          = 'g'; % we will draw one plume at a time, thus
                                   % specify a single color only
   
    opt = merge_options(opt, varargin{:});
    

    % plotting
    if isempty(h)
        h = figure;
    end
    num_tsteps = numel(plot_steps);

    
    for i = 1:numel(plot_steps)
        ix = plot_steps(i);
        subplot(1, num_tsteps, i);
        
        plume = reports(ix).sol.h;
        
        mapPlot(h, Gt, ...
                'plumes', plume, ...
                'plume_h_threshold', opt.plume_threshold, ...
                'plumeLineWidth', opt.plumeLineWidth, ...
                'plumeLineStyle', opt.plumeLineStyle, ...
                'casecolors', opt.plumeColor, ...
                'maplines', 0); % we don't want maplines, only plume outline
    end

end