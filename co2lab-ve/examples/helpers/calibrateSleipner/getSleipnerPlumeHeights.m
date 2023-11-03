function FH = getSleipnerPlumeHeights(varargin)

    opt = struct('year', 2001);
    opt = merge_options(opt, varargin{:});

    % We have 4 years worth of plume heights, taken from literature. The
    % heights are loaded from previously made .mat files, which are
    % available for download from the dataset manager: run "mrstDatasetGUI"
    % and click on "SleipnerPlumes".
    years   = [2001,2004,2006,2010];
    ind     = find(years == opt.year);
    if(isempty(ind))
        error('No plume for this year avilable');
    end
    datafiles = {   'ChadwickNoy_2001plume',...
                    'ChadwickNoy_2004plume',...
                    'ChadwickNoy_2006plume',...
                    'FurreEiken_2010plume'  };    
    datafile = fullfile(getDatasetPath('SleipnerPlumes'), [datafiles{ind},'.mat']);
    ok = exist(datafile,'file');
    if (ok>0)
        a   = load(datafile);
        FH  = @(x,y) interp2(a.ya, a.xa, a.H, y, x, 'linear', 0);
    else
        error('no data avilable');
    end
end