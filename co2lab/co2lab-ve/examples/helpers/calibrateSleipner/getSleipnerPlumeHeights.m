function FH = getSleipnerPlumeHeights(varargin)
%Undocumented Utility Function

%{
Copyright 2009-2024 SINTEF Digital, Mathematics & Cybernetics.

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
