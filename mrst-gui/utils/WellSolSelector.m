classdef WellSolSelector < UIItem
    properties
        Callback
        regCallback
        plotWellSolCallback
        nameSelector
        propSelector
        leftButton
        rightButton
        regionSwitch
        plotWSButton
        nameSubsetIx
        names
        props
        time
        name2prop
    end
    properties (Dependent)
        nameIx
        propIx
        panelNo
        curNames
        curProps
    end
    
    methods
        
        function s = WellSolSelector(namelist,proplist,time, varargin)
            opt = struct('Parent',          [], ...
                         'Callback',        @(src, event)disp('Callback unset'), ...
                         'regCallback',     @(src, event)disp('Region callback unset'), ...
                         'plotWellSolCallback',     @(src, event)disp('plotWellSol-callback unset'), ...
                         'Position',        [10 10 200 200], ...
                         'Visible',         'on', ...
                         'Title',           'Simulation output',...
                         'showPlotWSButton', true);
            
            [opt, extraOpt] = merge_options(opt, varargin{:});
            
            %[names, props, time, name2prop] = processSummary(smry);         
            wellHeader    = uicontrol('Parent', [], 'Style', 'text', ...
                'Value', [], 'String', 'Wells', 'Visible', 'off');
            propHeader   = uicontrol('Parent', [], 'Style', 'text', ...
                'Value', [], 'String', 'Properties', 'Visible', 'off');
            nameSelector = uicontrol('Parent', [], 'Style', 'listbox',  'Max', 2, 'Min', 0, ...
                                     'Value', [], 'String', namelist, 'Visible', 'off');
            propSelector = uicontrol('Parent', [], 'Style', 'listbox',  'Max', 2, 'Min', 0, ...
                                     'Value', [], 'String', proplist, 'Visible', 'off');  
                                 
            leftButton   = uicontrol('Parent', [], 'Style', 'radiobutton', 'String', 'Left panel', ...
                                     'Visible', 'off', 'Value', 1);
            rightButton  = uicontrol('Parent', [], 'Style', 'radiobutton', 'String', 'Right panel', ...
                                     'Visible', 'off');
            
            regionSwitch = uicontrol('Parent', [], 'Style', 'checkbox', 'String', 'Only wells for current region', ...
                                     'Visible', 'off'); 
             
            plotWSButton = uicontrol('Parent', [], 'Style', 'pushbutton', 'String', 'plotWellSols', ...
                                     'Visible', 'off');
            if opt.showPlotWSButton                   
                controls      = {{wellHeader,      propHeader},{nameSelector, propSelector},...
                    {leftButton, rightButton}, {regionSwitch}, {plotWSButton, []}};
                controlLayout = {[.5 .5],[.5, nan],[.5, nan], nan, [.5 nan]}; 
            else
                controls      = {{wellHeader,      propHeader},{nameSelector, propSelector},...
                    {leftButton, rightButton}, {regionSwitch}};
                controlLayout = {[.5 .5],[.5, nan],[.5, nan], nan}; 
            end
            
            
               
            s = s@UIItem('Parent', opt.Parent, 'controls', controls, ...
                         'controlWidths', controlLayout, 'Title', opt.Title, ...
                         'Position', opt.Position, 'Visible', 'off', extraOpt{:});
            
            [s.names, s.props] = deal(namelist, proplist);
            s.time = datetime(time, 'ConvertFrom', 'datenum');
            
            s.fixedHeight = true;
            s.nameSubsetIx = (1:numel(s.names))';
            % main callback
            s.Callback    = opt.Callback;
            s.regCallback = opt.regCallback;
            s.plotWellSolCallback = opt.plotWellSolCallback;
            % item callbacks
            nameSelector.Callback = @s.nameCallback;
            propSelector.Callback = @s.propCallback;
            [leftButton.Callback, rightButton.Callback] = deal(@s.buttonCallback);
            regionSwitch.Callback = @s.regionCallback;
            plotWSButton.Callback = @s.plotWellSolCallback;
            
            s.nameSelector = nameSelector;
            s.propSelector = propSelector;
            
            s.leftButton   = leftButton;
            s.rightButton  = rightButton;
            s.regionSwitch = regionSwitch;
            s.plotWSButton = plotWSButton;
            
            % set visible
            s.Visible = opt.Visible;            
        end
        
        function set.nameIx(s, val)
            s.nameSelector.Value = val;
        end
        function val = get.nameIx(s)
            val = s.nameSelector.Value;
        end
        
        function set.propIx(s, val)
            s.propSelector.Value = val;
        end
        function val = get.propIx(s)
            val = s.propSelector.Value;
        end
        
%         function set.panelNo(s, val)
%             if val == 1
%                 [s.leftButton.Value, s.rightButton.Value]  = deal(1, 0);
%             else
%                 [s.leftButton.Value, s.rightButton.Value]  = deal(0, 1);
%                 s.panelNo = 2;
%             end
%         end
        function val = get.panelNo(s)
            if s.leftButton.Value == 1
                val = 1;
            else
                val = 2;
            end
        end
        
        function set.curNames(s, val)
            ix = false(size(s.names));
            for k = 1:numel(val)
                ix = ix | strcmp(val{k}, s.names);
            end
            s.nameIx = find(ix);
        end
        function val = get.curNames(s)
            val = s.nameSelector.String(s.nameIx);
        end
        
        function set.curProps(s, val)
            ix = false(size(s.props));
            for k = 1:numel(val)
                ix = ix | strcmp(val{k}, s.props);
            end
            s.propIx = find(ix);
        end
        function val = get.curProps(s)
            val = s.propSelector.String(s.propIx);
        end        
           
        function nameCallback(s, src, event)
            %s.propSelector.Value  = [];
            %s.propSelector.String = s.props(any(s.name2prop(s.nameSubsetIx(s.nameIx),:), 1));
            s.Callback(src, event);
        end
        
        function propCallback(s, src, event)
            s.Callback(src, event);
        end
        
        function buttonCallback(s, src, event)
            if strncmp(src.String, 'Left', 4)   == 1
                s.rightButton.Value = abs(s.leftButton.Value-1);
            else
                s.leftButton.Value = abs(s.rightButton.Value-1);
            end
        end
        
        function regionCallback(s, src, event)
            s.nameSelector.Value  = [];
            s.propSelector.Value  = [];
            %s.propSelector.String = {};
            if s.regionSwitch.Value == 1
                s.regCallback(src, event)
            else
                s.nameSubsetIx = (1:numel(s.names))';
                s.nameSelector.String = s.names(s.nameSubsetIx);
            end
        end
        
        function setWellSubset(s, nms)
            [s.nameSelector.Value, s.propSelector.Value]  = deal([]);
            ix = false(numel(s.names), 1);
            for k = 1:numel(nms)
                nm = nms{k};
                ix = ix | strncmp(nm, s.names, 8);
            end
            s.nameSubsetIx = find(ix);
            s.nameSelector.String = s.names(s.nameSubsetIx);
        end
    end
end

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

% function [names, props, time, name2prop] = processSummary(smry)
% 
% if isfield(smry, 'MRST')
%     [names, props, time, name2prop] = processSummaryMRST(smry);
% else
%     [names, props, time, name2prop] = processSummaryECLIPSE(smry);
% end
% end


% function [names, props, time, name2prop] = processSummaryMRST(smry)
% 
%  % Use this for field name because any name could be used to get time and
%  % timestep data but this is unlikely to be the name of an actual field in 
%  % an MRST wellSol.
% specFld = ':+:+:+:+';
% time  = datenum(smry.STARTDAT) + smry.get(specFld, 'TIME', ':');
% names = smry.WGNAMES;
% props = smry.KEYWORDS(1:end-2);
% numwells = numel(smry.WGNAMES);
% numkws = numel(props);
% name2prop = ones(numwells,numkws);
% end
% 
% 
% 
% function [names, props, time, name2prop] = processSummaryECLIPSE(smry)
% specFld = ':+:+:+:+';
% time  = datenum(smry.STARTDAT) + smry.get(specFld, 'TIME', ':');
% nmIx  = ~strcmp(specFld, smry.WGNAMES);
% names = smry.WGNAMES(nmIx);
% props = smry.KEYWORDS;
% % name/prop compatibility matrix
% M = sparse(smry.nInx, smry.kInx, 1);
% M(~nmIx, :) = [];
% name2prop = logical(M);
% end
