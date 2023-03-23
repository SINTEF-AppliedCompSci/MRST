function data = getWeatherData(varargin)
%Get historical weather and climate data from MET Norway's databases. The
%function uses the Frost API, see https://frost.met.no/howto.html for
%details.

% (C) Øystein Klemetsdal (2022)

    opt = struct('userID'      , []                 , ...
                 'elements'    , {{'temperature'    , ...
                                   'precipitation'}}, ...
                 'startDate'    , '2020-06-01'      , ...
                 'endDate'      , '2021-05-31'      );
    opt = merge_options(opt, varargin{:});
    
    cmd = makeCommand();
    usr = getUser(opt);
    url = getURL();
    src = getSource();
    tmv = getTimeInterval(opt);
    ems = getElements(opt);
    trs = getTimeResolution();
    [fnm_str, fnm] = getFileName(src, tmv, ems, trs);
    
    if ~exist(fnm, 'file')
        command = [cmd, ' ', usr, ' ', '''', url, ...
                         src, '&', tmv, '&', ems, '&', trs, '''', ' ', fnm_str];
        system(command);
    end
    data = readJsonFile(fnm);
end

%-------------------------------------------------------------------------%
function str = makeCommand()
    str = 'curl -X GET';
end

%-------------------------------------------------------------------------%
function str = getUser(opt)
    userID = opt.userID;
    if isempty(userID)
        fn = fullfile(mrstPath('ruden-geothermal'), 'data', 'frost_id.txt');
        if ~exist(fn, 'file')
            error(['User ID for Frost API must be provided eihter as '     , ...
                   'optional input argument `userID`, or by making a file ', ...
                   '`data/frost_id.txt` with the user ID on a single line' ]);
        end
        userID = fileread(fn);
    end
    str = ['--user ', userID, ':'];
end

%-------------------------------------------------------------------------%
function str = getURL()
    str = 'https://frost.met.no/observations/v0.jsonld?';
end

%-------------------------------------------------------------------------%
function str = getSource()
    str = 'sources=SN18700%3A0';
end

%-------------------------------------------------------------------------%
function str = getTimeInterval(opt)
    str = strcat('referencetime=', opt.startDate, '%2F', opt.endDate);
end

%-------------------------------------------------------------------------%
function str = getElements(opt)
    elements = cellfun(@(element) translateElement(element), ...
                        opt.elements, 'UniformOutput', false);
    elements = cellfun(@(element) replace(element, ' ', '%20'), ...
                        elements, 'UniformOutput', false);
    elements = strjoin(elements, '%2C');
    str = strcat('elements=', elements);
end

%-------------------------------------------------------------------------%
function str = getTimeResolution()
    str = 'timeresolutions=PT1H';
end

%-------------------------------------------------------------------------%
function [str, fn] = getFileName(src, tmv, ems, trs)
    fn = fullfile(mrstPath('ruden-geothermal'), 'data');
    fn = fullfile(fn, [src, '_', tmv, '_', ems, '_', trs, '.json']);
    str = ['--output ', '"', fn, '"'];
end

%-------------------------------------------------------------------------%
function data = readJsonFile(fnm)
    % Read json file
    out  = jsondecode(fileread(fnm));
    data = out.data;
    % Get field names for output struct
    names = getFieldNames(data);
    % Allocate array for values and observation times
    nv     = numel(names);
    nd     = numel(data);
    values = nan(nd, nv);
    time   = cell(nd, 1);
    % Loop thorugh all observations and read data
    for i = 1:nd
        obs         = data(i).observations;
        values(i,:) = cellfun(@(obs) readObservation(obs), obs);
        time{i}     = data(i).referenceTime;
    end
    % Make output struct
    values = [num2cell(values, 1), {time}];
    data   = cell2struct(values, [names, 'time_stamp'], 2);
end
   
%-------------------------------------------------------------------------%
function names = getFieldNames(data)
    obs = data(1).observations;
    n = numel(obs);
    names = cell(1,n);
    for i = 1:n
        h = translateElement(obs{i}.elementId);
        if isfield(obs{i}, 'level')
            h = [h, '_', num2str(obs{i}.level.value)]; %#ok
        end
        names{i} = h;
    end
end

%-------------------------------------------------------------------------%
function val = readObservation(obs)
    val = obs.value;
    switch obs.unit
        case 'm'
            u = meter;
        case 'mm'
            u = milli*meter;
        otherwise
            u = 1;
    end
    val = val*u;
end

%-------------------------------------------------------------------------%
function out = translateElement(in)
    dict = { ...
                'air_temperature'               , 'temperature'  ;
                'sum(precipitation_amount PT1H)', 'precipitation';
           };
    subs = strcmpi(in, dict);
    if ~any(subs(:)), return; end
    [r,c] = find(subs);
    out   = dict{r, mod(c,2)+1};
end