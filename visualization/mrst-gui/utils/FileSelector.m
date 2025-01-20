classdef FileSelector < UIItem
    properties
        Callback
        upDirButton
        pathSelector
        pathSelectCutton
        fileNameEdit
        fileExtentionSelector
        fileListSelector
        applyButton
        recursiveSwitch        
        settings
        minFigSize = [nan, nan];
    end
    properties (Dependent)
        selectedFiles
    end
    
    methods
        
        function s = FileSelector(varargin)
            opt = struct('Parent',          [], ...
                         'Callback',        @(src, event)disp('Callback unset'), ...
                         'Position',        [10 10 700 500], ...
                         'Visible',         'on', ...
                         'Title',           'Select files', ...
                         'path',            '', ...
                         'file',            '', ...    
                         'fileExtensions',  {{'all'}}, ...
                         'settingsFile',    fullfile(mrstOutputDirectory(), 'file_selector_settings.mat'), ... 
                         'saveSettings',    true, ...
                         'standalone',      true);
            
            [opt, extraOpt] = merge_options(opt, varargin{:});
            
            %pathText = uicontrol('Parent', [], 'Style', 'text', 'Value', [], ... 
            %                     'String', 'Path:', 'Visible', 'off');
            upDirButton  = uicontrol('Parent', [], 'Style', 'pushbutton',  'Max', 1, 'Min', 0, ...
                                    'Value', [], 'String', '^', 'Visible', 'off');
            pathSelector = uicontrol('Parent', [], 'Style', 'popup', ...
                                    'Value', [], 'String', {}, 'Visible', 'off', 'Tooltip', 'Select path');
            pathSelectButton = uicontrol('Parent', [], 'Style', 'pushbutton',  'Max', 1, 'Min', 0, ...
                                   'Value', [], 'String', 'select dir', 'Visible', 'off');
            recursiveSwitch = uicontrol('Parent', [], 'Style', 'checkbox', 'String', 'Recursive search', ...
                                        'Visible', 'on', 'Value', 1);            
            %fileText = uicontrol('Parent', [], 'Style', 'text', 'Value', [], ... 
            %                     'String', 'File pattern:', 'Visible', 'off', 'HorizontalAlignment', 'right');
            fileNameEdit = uicontrol('Parent', [], 'Style', 'edit',  'Max', 1, 'Min', 0, ...
                                     'Value', [], 'String', '*', 'Visible', 'off', 'Tooltip', 'Select filename pattern');
            fileExtentionSelector = uicontrol('Parent', [], 'Style', 'popup',  'Max', 2, 'Min', 0, ...
                                              'Value', 1, 'String', opt.fileExtensions, 'Visible', 'off');
            fileListSelector = uicontrol('Parent', [], 'Style', 'listbox',  'Max', 2, 'Min', 0, ...
                                         'Value', [], 'String', {}, 'Visible', 'off'); 
            applyButton = uicontrol('Parent', [], 'Style', 'pushbutton',  'Max', 1, 'Min', 0, ...
                                    'Value', [], 'String', 'apply', 'Visible', 'off');

            controls      = {{upDirButton, pathSelector, pathSelectButton}, ...
                             {fileNameEdit, fileExtentionSelector, [], recursiveSwitch}, ...
                             {fileListSelector}, ...
                             {applyButton}};
            controlLayout = {[0.04, .84, .12], [.45, .2, .05, nan], nan, .2}; 
            
            fig = [];
            if ~isempty(opt.Parent)
                fig = opt.Parent;
            elseif opt.standalone
                fig = figure;
            end
            
            s = s@UIItem('Parent', fig, 'controls', controls, ...
                         'controlWidths', controlLayout, 'Title', opt.Title, ...
                         'Position', opt.Position, 'Visible', 'off', extraOpt{:});

            s.Callback    = opt.Callback;
            % item callbacks
            upDirButton.Callback = @s.upDir;
            pathSelector.Callback = @s.update;
            pathSelectButton.Callback = @s.pathSelectCallback;
            recursiveSwitch.Callback = @s.update;
            fileNameEdit.Callback = @s.update;
            fileExtentionSelector.Callback = @s.update;
            applyButton.Callback = @s.Callback;

            s.upDirButton = upDirButton;
            s.pathSelector = pathSelector;
            s.pathSelectCutton = pathSelectButton;
            s.recursiveSwitch = recursiveSwitch;
            s.fileNameEdit = fileNameEdit;
            s.fileExtentionSelector = fileExtentionSelector;
            s.fileListSelector = fileListSelector;
            s.applyButton = applyButton;

            s = update_settings(s, opt);
            pathSelector.String = s.settings.recentPaths;
            pathSelector.Value = 1;

            if opt.standalone
                dp = get(0, 'DefaultFigurePosition');
                s.Parent.SizeChangedFcn = @s.updateFigure;
                s.minFigSize = round([.8*dp(3), .5*dp(4)]);
                s.Parent.CloseRequestFcn = @s.closeFigure;
            end
            % set visible
            s.Visible = opt.Visible;     
            s.update();    
        end
        
        function updateFigure(s, varargin)
            mrg = [5, 5];
            fsz = s.Parent.Position(3:4);
            s.Position = [mrg, max(s.minFigSize, fsz)-2*mrg];
        end

        function list = get.selectedFiles(s)
            list = s.fileListSelector.String(s.fileListSelector.Value);
        end

        function closeFigure(s, src, event)
            if s.settings.save
                tmp = s.settings;
                save(s.settings.settingsFile, '-struct', 'tmp');
            end
            delete(s.Parent);
        end

        function fl = lookupFiles(s, maxNum, pth, depth)
            maxDepth = 2;
            if nargin < 4
                pth = s.pathSelector.String{s.pathSelector.Value};
                depth = 0;
            end
            extIx = s.fileExtentionSelector.Value;
            ext = s.fileExtentionSelector.String{extIx};
            if strcmp(ext, 'all')
                ext = '.*';
            end
            pattern = s.fileNameEdit.String;
            if isempty(pattern)
                pattern = '*';
            end
            df = dir(fullfile(pth, sprintf('%s%s', pattern, ext)));
            fl = applyFunction(@(d)fullfile(d.folder, d.name), df);
            if numel(fl) > maxNum
                fl = fl(1:maxNum);
                return;
            end
            if s.recursiveSwitch.Value == 1 && depth < maxDepth
                dd = dir(pth);
                dd = dd([dd.isdir]);
                for k = 1:numel(dd)
                    if dd(k).name(1) ~= '.'
                        fl = vertcat(fl, s.lookupFiles(maxNum-numel(fl), fullfile(pth, dd(k).name), depth + 1)); %#ok
                    end
                end
            end
        end

        function update(s, src, event)
            maxNumFiles = 100;
            curFiles = s.selectedFiles;
            curSelect = (1:numel(curFiles));
            newList = union(curFiles, s.lookupFiles(maxNumFiles-numel(curFiles)), 'stable');
            s.fileListSelector.String = newList(1:min(numel(newList), maxNumFiles));
            s.fileListSelector.Value = curSelect;
        end

        function pathSelectCallback(s, src, event)
            str = s.pathSelector.String{s.pathSelector.Value};
            msg = 'Select base folder';
            if ~isempty(str) || isfolder(str)
                pth = uigetdir(str, msg);
            else
                pth = uigetdir(msg);
            end
            s.addPath(pth);
        end

        function upDir(s, src, event)
            curPath = s.pathSelector.String{s.pathSelector.Value};
            s.addPath(fileparts(curPath));
        end

        function addPath(s, pth)
            if isa(pth, 'char')
                s.pathSelector.String = unique(vertcat({pth}, s.pathSelector.String), 'stable');
                s.settings.recentPaths = s.pathSelector.String;   
                s.pathSelector.Value = 1;
                s.update();
            end
        end
    end
end



function s = update_settings(s, opt)
inputpath = '';
if ~isempty(opt.path)
    tmp = what(opt.path);
    if ~isempty(tmp)
        inputpath = tmp.path;
    else
        warning('Could not find path: %s\n', opt.path);
    end
end
fn = opt.settingsFile;
ok = true;
if isfile(fn)
    s.settings = load(fn);
    ok = (isfield(s.settings, 'settingsFile') && ischar(s.settings.settingsFile)) && ...
         (isfield(s.settings, 'recentPaths') && iscell(s.settings.recentPaths));
end
if ~isfile(fn) || ~ok
    s.settings = struct('settingsFile', opt.settingsFile, ...
                        'recentPaths', {{}});
    if ~ok
        warning('Unexpected format of settings in %s. Not used.\n', fn);
    end
end
s.settings.save = opt.saveSettings;
if ~isempty(inputpath)
    s.settings.recentPaths = unique(vertcat({inputpath}, s.settings.recentPaths), 'stable');
end
if isempty(s.settings.recentPaths)
    s.settings.recentPaths = {pwd};
end
end


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
