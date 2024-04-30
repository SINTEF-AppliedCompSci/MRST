classdef SummarySelectorExtended < UIItem
    properties
        Callback
        nameSelector
        propSelector
        nameSearch
        propSearch
        typeSelect
        clearButton
        nameSubsetIx
        names
        props
        time
        nameMap
        propMap
        name2prop
        singlePropIx
    end
    properties (Dependent)
        nameIx
        propIx
        panelNo
        curNames
        curProps
    end
    
    methods
        
        function s = SummarySelectorExtended(smry, varargin)
            opt = struct('Parent',          [], ...
                         'Callback',        @(src, event)disp('Callback unset'), ...
                         'Position',        [10 10 200 200], ...
                         'Visible',         'on', ...
                         'Title',           'Simulation output', ...
                         'caseNames',       {applyFunction(@(n)sprintf('Case %d', n), (1:numel(smry)))} );
            
            [opt, extraOpt] = merge_options(opt, varargin{:});
            if ~isempty(smry)
                 [names, props, time, name2prop] = processSummary(smry);
            else
                error('');
            end
            popupStr = {'all', 'WELL', 'GROUP', 'FIELD',  'REGION', 'AQUIFER', 'BLOCK', 'MISC'};
            
            typeSelect = uicontrol('Parent', [], 'Style', 'popupmenu',  'Max', 1, 'Min', 0, ...
                                   'Value', 1, 'String', popupStr, 'Visible', 'off');
            clearButton = uicontrol('Parent', [], 'Style', 'pushbutton',  'Max', 1, 'Min', 0, ...
                                   'Value', [], 'String', 'clear all', 'Visible', 'off');
            
            nameSearch = uicontrol('Parent', [], 'Style', 'edit',  'Max', 1, 'Min', 0, ...
                                     'Value', [], 'String', '', 'Visible', 'off');
            propSearch = uicontrol('Parent', [], 'Style', 'edit',  'Max', 1, 'Min', 0, ...
                                     'Value', [], 'String', '', 'Visible', 'off');
            nameSelector = uicontrol('Parent', [], 'Style', 'listbox',  'Max', 2, 'Min', 0, ...
                                     'Value', [], 'String', names, 'Visible', 'off');
            propSelector = uicontrol('Parent', [], 'Style', 'listbox',  'Max', 2, 'Min', 0, ...
                                     'Value', [], 'String', props, 'Visible', 'off');                     
                               
            controls      = {{typeSelect, clearButton}, {nameSearch, propSearch}, {nameSelector, propSelector}};
            controlLayout = {[.5, .5], [.5, .5], [.5, .5]}; 
            
            
            s = s@UIItem('Parent', opt.Parent, 'controls', controls, ...
                         'controlWidths', controlLayout, 'Title', opt.Title, ...
                         'Position', opt.Position, 'Visible', 'off', extraOpt{:});
            [s.names, s.props] = deal(names, props);
            [s.time, s.name2prop] = deal(time, name2prop);
            
            [nnm, nprp] = deal(numel(names), numel(props));
            s.nameMap = struct('sub',(1:nnm)' ,  'ix', ones(nnm,1));
            s.propMap = struct('sub',(1:nprp)', 'ix',  ones(nprp,1));
            s.singlePropIx = cellfun(@any, regexp(s.props, '^[^ACGLRSW]|^STEP'));
            
            s.fixedHeight = true;
            s.nameSubsetIx = (1:numel(s.names))';
            % main callback
            if ~isempty(smry)
                s.Callback    = opt.Callback;
                % item callbacks
                nameSelector.Callback = @s.nameCallback;
                propSelector.Callback = @s.propCallback;
                nameSearch.Callback   = @s.searchCallback;
                propSearch.Callback   = @s.searchCallback;
                typeSelect.Callback   = @s.typeCallback;
                clearButton.Callback  = @s.clearCallback;
                
                s.nameSelector = nameSelector;
                s.propSelector = propSelector;
                s.nameSearch   = nameSearch;
                s.propSearch   = propSearch;
                s.typeSelect   = typeSelect;
                s.clearButton  = clearButton;
            end
            
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
        
        function set.panelNo(s, val)
            if val == 1
                [s.leftButton.Value, s.rightButton.Value]  = deal(1, 0);
            else
                [s.leftButton.Value, s.rightButton.Value]  = deal(0, 1);
                s.panelNo = 2;
            end
        end
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
            s.nameMap.ix(s.nameMap.ix == 2) = 1;
            s.nameMap.ix(s.nameMap.sub(s.nameSelector.Value)) = 2;
            if isempty(s.propSearch.String) && ~isempty(s.nameSelector.Value)
                gix = any(s.name2prop(s.nameMap.ix == 2,:)', 2);
                s.setPropListSubset(gix);
            end
            s.Callback(src, event);
        end
        
        function propCallback(s, src, event)
            s.propMap.ix(s.propMap.ix == 2) = 1;
            s.propMap.ix(s.propMap.sub(s.propSelector.Value)) = 2;
            % if there's a single name-match, select name auto
            isSingleSelected = (s.propMap.ix == 2) & s.singlePropIx;
            if any(isSingleSelected)
                gix = any(s.name2prop(:,isSingleSelected), 2);
                s.setNameListSelect(gix);
            end
            s.Callback(src, event);
        end
        
        function searchCallback(s, src, event)
            namePat = s.nameSearch.String;
            propPat = s.propSearch.String;
            if isempty(namePat)
                nix = true(size(s.names));
            else
                if isvarname(namePat)
                    namePat = ['^', namePat];
                end
                nix = cellfun(@any, regexp(s.names, namePat));
            end
            if isempty(propPat)
                pix = true(size(s.props));
            else
                % reset type-selection
                s.typeSelect.Value = 1;
                if isvarname(propPat)
                    propPat = ['^', propPat];
                end
                pix = cellfun(@any, regexp(s.props, propPat));
            end
            pix2 = any(s.name2prop(nix, :)', 2);% & s.propMap.ix > 0;
            nix2 = any(s.name2prop(:, pix), 2);%  & s.nameMap.ix > 0;
            s.setNameListSubset(nix & nix2);
            s.setPropListSubset(pix & pix2);
            s.Callback(src, event); 
        end
        
        function searchForNames(s, src, event)
            pat  = s.nameSearch.String;
            % if no special char, assume "starts with"
            if isvarname(pat)
                pat = ['^', pat];
            end
            if ~isempty(pat)
                mix = cellfun(@any, regexp(s.names, pat));
            else
                mix = true(size(s.names));
            end
            s.setNameListSubset(mix);
            pix = s.propMap.ix > 0;
            if true%~isempty(s.propSearch.String)
                gix = double(any(s.name2prop(mix, :), 1));
                pix = pix & gix(:);
            end
            s.setPropListSubset(pix);
            s.Callback(src, event);
        end
        
        function searchForProps(s, src, event)
            pat  = s.propSearch.String;
            if isvarname(pat)
                pat = ['^', pat];
            end
            if ~isempty(pat)
                mix = cellfun(@any, regexp(s.props, pat));
            else
                mix = true(size(s.props));
            end
            s.setPropListSubset(mix);
            nmix = s.nameMap.ix > 0;
            if true%~isempty(s.nameSearch.String)
                gix = double(any(s.name2prop(:, mix)', 1));
                nmix = nmix & gix(:);
            end
            s.setNameListSubset(nmix);
            s.Callback(src, event);
        end
        
        function s = setPropListSubset(s, gix)
            globSel = s.propMap.ix == 2;
            gix = double(gix);
            gix(globSel) = 2;
            s.propMap.ix  = gix;
            s.propMap.sub = find(gix);
            s.propSelector.String = s.props(s.propMap.sub);
            s.propSelector.Value = find(s.propMap.ix(s.propMap.sub) == 2);
        end
       
        
        function s = setNameListSubset(s, gix)
            globSel = s.nameMap.ix == 2;
            gix = double(gix);
            gix(globSel) = 2;
            s.nameMap.ix  = gix;
            s.nameMap.sub = find(gix);
            s.nameSelector.String = s.names(s.nameMap.sub);
            s.nameSelector.Value = find(s.nameMap.ix(s.nameMap.sub) == 2);
        end
        
        function s = clearCallback(s, src, event)
            s.typeSelect.Value = 1;
            s.nameSearch.String = '';
            s.propSearch.String = '';
            s.nameSelector.Value = [];
            s.nameCallback(src, event);
            s.propSelector.Value = [];
            s.propCallback(src, event);
            s.searchCallback(src, event);
            
        end
        
        function s = typeCallback(s, src, event)
            tp = s.typeSelect.String{s.typeSelect.Value};
            str = '';
            switch tp
                case 'WELL'
                    str = '^[WCS]';
                case 'GROUP'
                    str = '^G';
                case 'FIELD'
                    str = '^F';
                case 'REGION'
                    str = '^R';
                case 'AQUIFER'
                    str = '^A';                
                case 'MISC'
                    str = '^[^ABCFGLRSW]|^STEP';
                case 'BLOCK'
                    str = '^B';                  
            end
            s.propSearch.String = str;
            s.searchCallback(src, event);
        end
        
        function s = setNameListSelect(s, gix)
            s.nameMap.ix(gix) = 2;
            s.nameMap.sub = find(s.nameMap.ix);
            s.nameSelector.String = s.names(s.nameMap.sub);
            s.nameSelector.Value = find(s.nameMap.ix(s.nameMap.sub) == 2);
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
            s.propSelector.String = {};
            if s.regionSwitch.Value == 1
                s.regCallback(src, event)
            else
                s.nameSubsetIx = (1:numel(s.names))';
                s.nameSelector.String = s.names(s.nameSubsetIx);
            end
        end
    end
end

function [names, props, time, name2prop] = processSummary(smry)
names  = applyFunction(@(s)s.WGNAMES, smry);
nInx   = applyFunction(@(s)s.nInx, smry);
[names, nInx] = uniqueStringIndex(names, nInx);
props  = applyFunction(@(s)s.KEYWORDS, smry);
kInx   = applyFunction(@(s)s.kInx, smry);
[props, kInx] = uniqueStringIndex(props, kInx);
M = sparse(nInx, kInx, 1);
name2prop = logical(M);
time = cell(1, numel(smry));
for k  = 1:numel(smry)
    specFld = smry{k}.getNms('TIME');
    tmp  = datenum(smry{k}.STARTDAT) + smry{k}.get(specFld, 'TIME', ':');
    time{k} = datetime(tmp, 'ConvertFrom', 'datenum');
end
end

function [ss, ii] = uniqueStringIndex(ss, ii)
cumn = cumsum( [1, cellfun(@numel, ss)]);
[ss, ~, nix] = unique(vertcat(ss{:}));
for k = 1:numel(ii)
    tmp = nix(cumn(k):(cumn(k+1)-1));
    ii{k} = tmp(ii{k});
end
ii = vertcat(ii{:});
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
