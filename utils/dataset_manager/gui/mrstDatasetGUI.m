function varargout = mrstDatasetGUI()
%Open dataset management user interface
%
% SYNOPSIS:
%   h = mrstDatasetGUI()
%
% DESCRIPTION:
%   The dataset GUI is a user interface, allowing the user to see which
%   datasets are known to MRST, visit their webpages, download and manage
%   them. The function is built upon the dataset library.
%
% RETURNS:
%   h   - Handle to panel figure.
%
%
% SEE ALSO:
%   `getDatasetInfo`, `getDatasetPath`, `listDatasetExamples`

%{
Copyright 2009-2021 SINTEF Digital, Mathematics & Cybernetics.

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

    persistent h
    names = {};

    if ~isempty(h) && ishandle(h)
        % Avoid creating more than one window, since we can't rely on close all
        % to get rid of them.
        close(h)
    end

    figNo = 100;
    while ishandle(figNo)
        figNo = figNo + 1;
    end
    h = figure(figNo);
    set(h,     'Toolbar',       'none', ...
               'NumberTitle',   'off', ...
               'MenuBar',       'none');
    [info, present, sh] = deal([]);
    function startupGUI()
        [info, present] = getAvailableDatasets();
        names = arrayfun(@(x) x.name, info, 'UniformOutput', false);
        if ~isempty(sh)
            drawSelection(sh, [])
        end
    end
    startupGUI();

    sh = uicontrol(h, 'Style', 'listbox', 'units', 'normalized', ...
                 'Callback', @drawSelection, ...
                 'String', names,  'position', [0, 0, .3, 1]);
    ph = uipanel(h, 'position',[0.3, 0, .7, 1], 'backgroundcolor', 'w');
    
    function drawSelection(src, event)
        delete(get(ph, 'Children'));
        
        current = get(src, 'Value');
        I = info(current);
        
        set(h, 'Name', ['MRST Dataset manager: ', names{current}]);
        
        % Check for image. If we have a image, we allocate some room and
        % place it in the panel.
        hasImage = ~isempty(I.image);
        if hasImage
            voffset = 0;
            ax = axes('parent', ph, 'units', 'normalized', ...
                        'position', [0, .5, 1, .5]);
            c = imread(I.image);
            image(c, 'parent', ax);
            axis equal tight off
        else
            voffset = .5;
        end
        
        txt = [I.modelType, ' model'];
        if ~isempty(I.cells)
            txt = [txt, ' with ', num2str(I.cells) ' cells'];
        end
        % Use cell array for line breaks in edit box
        txt = {txt, getContents(I), I.description};
        if ~isempty(I.note)
           txt = [ txt, {''}, reshape(I.note, 1, []) ];
        end

        uicontrol(ph, 'units', 'normalized', ...
                       'horizontalAlignment', 'left', ...
                       'Max', 1000, 'Min', 0, ...
                       'position', [0, .2, 1, .30 + voffset],...
                       'BackgroundColor', 'w', ...
                       'Style', 'Edit', 'String', txt);
        
        canDownload = datasetHasCustomDownloadFcn(I) || ...
                      datasetHasValidFileURL     (I);

        onDisk = present(current);
        
        if onDisk
            statusstr = 'Present';
            fc = [0 1 0];
        else
            if canDownload
                statusstr = 'Available for download';
                fc = [1, 165/255, 0];
            elseif isempty(I.website);
                statusstr = 'Not available for download';
                fc = [1, 0, 0];
            else
                statusstr = 'See website';
                fc = [1, 0, 0];
            end
        end
        uicontrol(ph, 'units', 'normalized', ...
                       'horizontalAlignment', 'center', ...
                       'position', [0, .1, 1, .08],...
                       'ForegroundColor', fc, ...
                       'backgroundcolor', 'w', ...
                       'FontSize', 18, ...
                       'Style', 'Text', 'String', statusstr);
        
        st = {'off', 'on'};
        addbutton = @(x, label, active, varargin) uicontrol(ph,...
                                   'Units', 'normalized', ...
                                   'Enable', st{1 + active}, ...
                                   'Position', [x, 0, 0.2, .1], ...
                                   'Style', 'pushbutton', ...
                                   'String', label, ...
                                   varargin{:});
        
        addbutton(.0, 'Download', canDownload,...
                      'Callback', @(src, event) getDataset(I, present(current)))
        addbutton(.2, 'Delete', onDisk, 'Callback', @(src, event) deleteDataset(I))
        addbutton(.4, 'Webpage', ~isempty(I.website), 'Callback', @(src, event) web(I.website))
        addbutton(.6, 'List files', present(current), 'Callback', @(src, event) listfiles(I))
        addbutton(.8, 'Examples', ~isempty(I.examples), 'Callback', @(src, event) listDatasetExamples(I))
    end

    drawSelection(sh, [])
    
    function deleteDataset(info)
        pth = fullfile(mrstDataDirectory(), info.name);
        status = questdlg(['Please confirm deletion of dataset. This will delete all files under', ...
                           pth, ...
                           ', including any user modifications!'], ...
                           'title','Yes','No','No');
        switch(lower(status))
            case 'yes'
                rmdir(pth, 's');
                startupGUI()
            case 'no'
                % Do nothing
        end
    end

    function getDataset(info, present)
        if present
            helpdlg([info.name, ' dataset already downloaded!']);
            return
        end

        status = questdlg(['Would you like to download the ', info.name, ...
                           ' dataset from the Internet? Download size: ', num2str(round(info.filesize)) ' MB.'],...
                           'title','Yes','No','Yes');
        switch(lower(status))
            case 'yes'
                downloadDataset(info.name, false);
                msgbox('Successfully downloaded dataset!');
                startupGUI()
            case 'no'
                % Do nothing
        end
    end
   if nargout > 0
       varargout{1} = h;
   end
   set(h, 'HandleVisibility',  'Callback')
end

function listfiles(info)
    pth = fullfile(mrstDataDirectory(), info.name);
    disp('**** Dataset installed in :');
    disp(['**** ', pth]);
    disp('Directory listing:');
    ls(pth)
end


function included = getContents(I)
    stuff = [I.hasGrid, I.hasRock, I.hasFluid];
    stuffn = {'Grid', 'Rock', 'Fluid'};
    nn = nnz(stuff);

    count = 0;
    included = '';
    for i = 1:numel(stuff)
        if ~stuff(i)
            continue
        end
        if count > 0
            included = [included, lower(stuffn{i})];
        else
            included = ['Model contains: ', stuffn{i}];
        end

        count = count + 1;

        if count == nn
            included = [included, '.'];
        elseif count == nn - 1
            included = [included, ' and '];
        else
            included = [included, ', '];
        end
    end

end
