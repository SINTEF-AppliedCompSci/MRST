function [wells, states] = readJutulOutput(pth, varargin)
    opt = struct('wells', true, 'states', true, 'readStates', true);
    opt = merge_options(opt, varargin{:});
    [filepath, name, ext] = fileparts(pth);

    fldr = sprintf('%s_output_mrst', name);
    if opt.states
        states = readStates(filepath, fldr, opt.readStates);
    else
        states = [];
    end
    well_path = fullfile(filepath, fldr, 'wells.mat');
    if exist(well_path, 'file')
        wells = load(well_path);
    else
        wells = [];
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