function output = plotSimplePDFCDF(x, varargin)
% Simple plot of Probability Distribution Function (PDF) and
% Cumulative Distribution Function (CDF)

    opt.numBins = [];
    opt.xlabel = 'CO2 (Mt)'; % default
    opt.titleOn = true;
    opt.textboxOn = false;
    opt = merge_options(opt, varargin{:});

    figure, hold on

    yyaxis left
    if isempty(opt.numBins)
        histogram(x)
    else
        histogram(x,opt.numBins) 
    end
    xlabel(opt.xlabel)
    ylabel('Frequency')

    yyaxis right
    plot(sort(x), (1:numel(x))./numel(x), 'LineWidth',2) % y-axis goes from 0 to 1
    ylabel('Probability')

    meanX = mean(x);
    standardDev = std(x); % normalizes by N-1
    variance = var(x); % normalizes by N-1

    sorted_x = sort(x);
    r = numel(x);
    P10inx = round(r*0.1); if P10inx == 0, P10inx = 1; end;
    P50inx = round(r*0.5);
    P90inx = round(r*0.9);
    P10val = sorted_x(P10inx);
    P50val = sorted_x(P50inx);
    P90val = sorted_x(P90inx);

    if opt.textboxOn
        ha = annotation('textbox',[.15 .6 .175 .3], ...
            'String',  {['mean=',num2str(meanX,3)]; ...
                        ['SD=',num2str(standardDev,3)]; ...
                        ['var=',num2str(variance,3)]} );
        ha.EdgeColor = 'none';
    end
    % Or include in title:
    if opt.titleOn
        title({['mean=',num2str(meanX,3),', SD=',num2str(standardDev,3),', var=',num2str(variance,3)]; ...
           ['P10=',num2str(P10val,3),', P50=',num2str(P50val,3),', P90=',num2str(P90val,3)]}   );
    end
    box on;

    
    % Possible output:
    output.mean = mean(x);
    output.median = median(x);
    output.oneStd = std(x);
    output.cov = std(x) / mean(x) * 100;
    output.P10 = P10val;
    output.P50 = P50val;
    output.P90 = P90val;
    output.min = min(x);
    output.max = max(x);
    output.spread = max(x) - min(x);


end
