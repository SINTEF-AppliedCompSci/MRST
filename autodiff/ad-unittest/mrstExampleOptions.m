function opt = mrstExampleOptions(filename)
%Parse MRST_TEST_OPTIONS annotation block from an example file.
%
% SYNOPSIS:
%   opt = mrstExampleOptions(filename)
%
% DESCRIPTION:
%   Reads the structured comment block ``MRST_TEST_OPTIONS`` near the top of
%   a MATLAB example file and returns the parsed key/value pairs as a struct.
%
%   The annotation block has the following form (inside a ``%{ ... %}``
%   block)::
%
%       %{
%       MRST_TEST_OPTIONS
%       timeout:     120        % seconds; 0 = no limit
%       interactive: false      % true = skip in automated runs
%       tags:        slow, data-required   % comma-separated tags
%       %}
%
%   All fields are optional.  Defaults are:
%     timeout     - 0 (no limit)
%     interactive - false
%     tags        - {} (no tags)
%
% PARAMETERS:
%   filename - Full path to a MATLAB script (.m file), or just the name
%              without extension.  If no path is given the file is located
%              via the MATLAB path.
%
% RETURNS:
%   opt - Struct with fields:
%           timeout     (double, seconds)
%           interactive (logical)
%           tags        (cell array of char)
%
% EXAMPLES:
%   opt = mrstExampleOptions('myExampleScript');
%   if opt.interactive
%       disp('This example requires user interaction');
%   end
%
% SEE ALSO:
%   MRSTExampleTests, getExampleIntegrationTestSuiteMRST

%{
Copyright 2009-2026 SINTEF Digital, Mathematics & Cybernetics.

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

    % --- Default options -----------------------------------------------
    opt = struct('timeout'    , 0    , ...
                 'interactive', false, ...
                 'tags'       , {{}} );

    % --- Resolve filename -----------------------------------------------
    if ~ischar(filename)
        return
    end
    % Strip .m extension if present
    [~, ~, ext] = fileparts(filename);
    if strcmpi(ext, '.m')
        filename = filename(1:end-2);
    end
    % Find on MATLAB path if not an absolute path
    if ~isAbsolutePath(filename)
        fullp = which([filename, '.m']);
        if isempty(fullp)
            return
        end
        filename = fullp(1:end-2);
    end
    filepath = [filename, '.m'];
    if ~exist(filepath, 'file')
        return
    end

    % --- Read file and search for annotation block ----------------------
    fid = fopen(filepath, 'r');
    if fid < 0
        return
    end
    lines = {};
    while ~feof(fid)
        lines{end+1} = fgetl(fid); %#ok<AGROW>
    end
    fclose(fid);

    % Search within the first 50 lines for the MRST_TEST_OPTIONS block
    limit = min(50, numel(lines));
    inBlock  = false;
    blockOpen = false;   % inside a %{ ... %} comment
    kvLines  = {};

    for i = 1:limit
        ln = strtrim(lines{i});
        if ~inBlock
            if strcmp(ln, '%{')
                blockOpen = true;
            elseif blockOpen && strcmpi(ln, 'mrst_test_options')
                inBlock = true;
            elseif blockOpen && strcmp(ln, '%}')
                blockOpen = false;
            end
        else
            if strcmp(ln, '%}')
                break
            end
            kvLines{end+1} = ln; %#ok<AGROW>
        end
    end

    if isempty(kvLines)
        return
    end

    % --- Parse key: value lines -----------------------------------------
    for i = 1:numel(kvLines)
        ln = kvLines{i};
        % Strip trailing inline comment (% ...)
        cIdx = strfind(ln, '%');
        if ~isempty(cIdx)
            ln = strtrim(ln(1:cIdx(1)-1));
        end
        if isempty(ln), continue; end

        colonIdx = strfind(ln, ':');
        if isempty(colonIdx), continue; end

        key   = lower(strtrim(ln(1:colonIdx(1)-1)));
        value = strtrim(ln(colonIdx(1)+1:end));

        switch key
            case 'timeout'
                v = str2double(value);
                if ~isnan(v)
                    opt.timeout = v;
                end
            case 'interactive'
                opt.interactive = any(strcmpi(value, {'true','1','yes','on'}));
            case 'tags'
                parts = strsplit(value, ',');
                opt.tags = cellfun(@strtrim, parts, 'UniformOutput', false);
                opt.tags = opt.tags(~cellfun(@isempty, opt.tags));
        end
    end
end

% =========================================================================

function tf = isAbsolutePath(p)
    if ispc
        tf = numel(p) >= 3 && p(2) == ':' && (p(3) == '\' || p(3) == '/');
    else
        tf = ~isempty(p) && p(1) == '/';
    end
end
