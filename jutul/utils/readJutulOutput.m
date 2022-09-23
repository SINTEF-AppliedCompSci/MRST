function [wells, states] = readJutulOutput(pth, varargin)
% Read Jutul output from a path given by writeJutulInput once simulated.

%{
Copyright 2009-2022 SINTEF Digital, Mathematics & Cybernetics.

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

    opt = struct('wells', true, 'states', true, ...
                 'readStates', true, 'wellSol', true, ...
                 'error', false);
    opt = merge_options(opt, varargin{:});
    [filepath, name, ext] = fileparts(pth);

    [states, wells] = deal([]);
    ok = false;
    fldr = sprintf('%s_output_mrst', name);
    if opt.states
        states = readStates(filepath, fldr, opt.readStates);
        ok = ~isempty(states);
    end
    well_path = fullfile(filepath, fldr, 'wells.mat');
    if opt.wells
        if exist(well_path, 'file')
            wells = load(well_path);
            if opt.wellSol
                wells = convertJutulWellSols(wells);
                if opt.states
                    for i = 1:numel(states)
                        % states{i}.wellSol = wells{i};
                    end
                end
            end
        else
            ok = false;
        end
    end
    if ~ok
        msg = 'No Jutul output found. Case may not be simulated?';
        if opt.error
            error(msg);
        else
            warning(msg);
        end
    end
end

function states = readStates(directory, fldr, readStates)
    h = ResultHandler('writeToDisk', true, 'dataDirectory', directory, 'dataFolder', fldr);
    ns = h.numelData();
    if readStates
        % Deserialize all the states
        states = cell(ns, 1);
        for i = 1:ns
            states{i} = h{i};
        end
    else
        % We just want the handler
        states = h;
    end
end
