function [] = setUIStyle(d, styleName, layoutOpt)
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

if nargin < 3
    layoutOpt = {};
end
if nargin < 2
    styleName = 'default';
end
% definne some colors
standard = [.94 .94 .94];
blue     = [0.1328 0.2109 0.3496];

% defaults
props = struct('BackgroundColor',        standard, ...
               'ForegroundColor',            blue, ...
               'titleColor',                 blue, ...
               'dummyColor',             standard, ...
               'FontName',            'Helvetica', ...
               'titleFontName',       'Helvetica', ...
               'FontSize',                     10, ...
               'titleFontSize',                10, ...
               'FontWeight',             'normal', ...
               'titleFontWeight',        'normal', ...
               'titleFontAngle',         'normal', ...
               'FontAngle',              'normal', ...         
               'HorizontalAlignment',      'left');

panelProps = struct('HighlightColor',       blue, ...
                    'ShadowColor',      standard, ...
                    'BorderType',         'none', ...
                    'BorderWidth',            0, ...
                    'TitlePosition',   'lefttop');
                
dummyPanelProps = panelProps;                
   
% set default parameters
if isa(d, 'UIItem')
    d.layout.params = getDefaultsItemParams();
elseif isa(d, 'UIMenu')
    d.layout.params = getDefaultsMenuParams();
end

% reset props according to style
switch styleName
    % DEFAULT (single color, same background color) -------------------------------------------------------------
    case 'default'
        % all is allready set, just make sure there is no dummy-panel
        delete(d.dummyPanel)
        d.dummyPanel = [];
    % CLEAN (single color, menu extent shown as line to the left) ----------------------------------------------- 
    case 'clean'
        if ~isa(d.dummyPanel, 'matlab.ui.container.Panel') || ~isvalid(d.dummyPanel)
            d.dummyPanel = uipanel('Units', 'pixels', 'Position', d.panel.Position, 'Visible', 'on');
        end
        d.dummyPanel.Parent = d.panel.Parent;
        % position behind
        uistack(d.dummyPanel, 'bottom');
        dummyPanelProps.BorderType   = 'line';
        dummyPanelProps.BorderWidth  = 1;
        dummyPanelProps.Title =  '    .';
        d.layout.params.margins(1)      =  3;
        d.layout.params.outerMargins(1) =  3;
    %BOARDER (single color, standard boarder) ------------------------------------------------------------------ 
    case 'boarder'
        delete(d.dummyPanel)
        d.dummyPanel = [];
        panelProps.BorderType   = 'line';
        panelProps.BorderWidth  = 1;
    % BOX (single color, box around each item) ------------------------------------------------------------------ 
    case 'box'
        if ~isa(d.dummyPanel, 'matlab.ui.container.Panel') || ~isvalid(d.dummyPanel)
            d.dummyPanel = uipanel('Units', 'pixels', 'Position', d.panel.Position, 'Visible', 'on');
        end
        d.dummyPanel.Parent = d.panel.Parent;
        % position behind
        uistack(d.dummyPanel, 'bottom');
        dummyPanelProps.BackgroundColor = blue;
        
        d.layout.params.outerMargins    = [1 0 1 1];
        if isa(d, 'UIMenu') && d.level == 1
            d.layout.params.outerMargins    = [1 1 1 1];
        end
        if isa(d, 'UIMenu')
            d.layout.params.vskip = 5;
            d.layout.params.margins    = [4 0 0 0];
        else
            d.layout.params.margins    = [4 4 5 0];
        end
    %INVERT (single color, switch colors of title/title-background) --------------------------------------------- 
    case 'invert'
        if ~isa(d.dummyPanel, 'matlab.ui.container.Panel') || ~isvalid(d.dummyPanel)
            d.dummyPanel = uipanel('Units', 'pixels', 'Position', d.panel.Position, 'Visible', 'on');
        end
        d.dummyPanel.Parent = d.panel;
        uistack(d.dummyPanel, 'bottom');
        %uistack(d.dummyPanel, 'up');
        d.dummyStack = 'top';
        if isa(d, 'UIItem')
            props.titleBackgroundColor = blue;
            props.BackgroundColor = standard;
            props.titleColor = standard;
        else
            props.BackgroundColor = blue;
            props.ForegroundColor = standard;
        end
        if isa(d, 'UIMenu') && d.level == 1
            d.layout.params.dummyOuterMargins    = [2 2 2 2];
        else
            d.layout.params.dummyOuterMargins    = [2 2 2 2];
        end
        if isa(d, 'UIMenu')
            d.layout.params.margins              = [0 0 0 0];
        else
            d.layout.params.margins              = [5 5 1 1];
        end
    %DEBUG (layout for debug purposes) ------------------------------------------------------------------ 
    case 'debug'
        if ~isa(d.dummyPanel, 'matlab.ui.container.Panel') || ~isvalid(d.dummyPanel)
            d.dummyPanel = uipanel('Units', 'pixels', 'Position', d.panel.Position, 'Visible', 'on');
        end
        d.dummyPanel.Parent = d.panel.Parent;
        % position behind
        uistack(d.dummyPanel, 'bottom');
        dummyPanelProps.BackgroundColor = blue;
        d.layout.params.outerMargins    = [1 1 1 1];
        if isa(d, 'UIMenu')
            d.layout.params.margins         = [1 1 1 1];
        else
            d.layout.params.margins         = [1 1 1 1];
        end
    otherwise
        warning('Unkown style: %s', styleName);  
end

% merge in user-given properties and set props
[props, ~] = merge_options(props, layoutOpt{:});
setProps(d, props);
setProps(d.panel, panelProps);
setProps(d.dummyPanel, dummyPanelProps)

% update layout-field and set titleHeight
d.layout.styleName = styleName;
d.layout.opt  = layoutOpt; 
yMrg = sum(d.layout.params.outerMargins(3:4));
% Measure actual panel title bar height dynamically so the layout adapts to
% the MATLAB version (the title bar pixel height changed between R2024b and
% R2025a).  Fall back to the font-size formula when the measurement is not
% yet available (e.g. panel not yet drawn).
d.titleHeight = getPanelTitleHeight(d.panel) + yMrg;
end
 
function setProps(d, props)
pnms = fieldnames(props);
for k = 1:numel(pnms)
    if isprop(d, pnms{k})
        d.(pnms{k}) = props.(pnms{k});
    end
end
end

function params = getDefaultsItemParams()
    params = struct('margins',            [5 5 2 2], ...  % left, right, bottom, top
                    'outerMargins',       [0 0 0 0], ...  % left, right, bottom, top 
                    'dummyOuterMargins',  [0 0 0 0], ...  % 
                    'vskip',                      5, ...  % vertical space between items
                    'hskip',                      2, ...  % horizontal space between items
                    'lineHeightAt10',            20, ...  % linehight at fonsize 10 pt
                    'editWidth',                 60);     % fixed (pt)
end

function params = getDefaultsMenuParams()
    params = struct('margins',            [5 5 0 0], ...  % left, right, bottom, top
                    'outerMargins',       [0 0 0 0], ...  % left, right, bottom, top 
                    'dummyOuterMargins',  [0 0 0 0], ...  %
                    'vskipMenu',                  1, ...
                    'vskipItem',                  5);       
end

function h = getPanelTitleHeight(panel)
% Reliably measure the uipanel title bar height across MATLAB versions.
%
% Strategy (in order of preference):
%  1. Create a throw-away off-screen test panel that has a visible border
%     (BorderType='line') so InnerPosition accurately excludes the title bar.
%     Force TitlePosition='lefttop' so the title is at the top regardless of
%     the MATLAB version default.  Use drawnow('nocallbacks') to flush the
%     layout engine before reading InnerPosition.
%  2. Version-based fallback: R2025a+ uses a taller title bar (~24 px at 10
%     pt) than R2024b and earlier (~18 px at 10 pt).
    try
        % Create an off-screen test panel so the measurement causes no flicker.
        % BorderType='line' gives a predictable InnerPosition:
        %   InnerPosition.y   = borderWidth  (bottom border)
        %   InnerPosition.h   = panelH - titleBarH - borderWidth  (bottom border only,
        %                       title replaces the top border)
        % => titleBarH = panelH - InnerPosition.y - InnerPosition.h
        tmpPanel = uipanel('Parent',        panel.Parent, ...
                           'Units',         'pixels', ...
                           'FontName',      panel.FontName, ...
                           'FontSize',      panel.FontSize, ...
                           'TitlePosition', 'lefttop', ...
                           'BorderType',    'line', ...
                           'BorderWidth',   1, ...
                           'Title',         'T', ...
                           'Position',      [0 -9999 200 200]);
        drawnow('nocallbacks');
        ip = tmpPanel.InnerPosition;
        pp = tmpPanel.Position;
        delete(tmpPanel);
        h = pp(4) - ip(2) - ip(4);
        if h > 0 && h <= 60
            return;
        end
    catch
        % ignore and fall through to fallback
    end
    % Version-based fallback.  R2025a changed the default UI font metrics,
    % making the panel title bar taller.
    try
        v    = version('-release');           % e.g. '2024b' or '2025a'
        year = str2double(v(1:4));            % numeric year
    catch
        year = 0;
    end
    if year >= 2025
        h = ceil(24 * panel.FontSize / 10);  % R2025a and later
    else
        h = ceil(18 * panel.FontSize / 10);  % R2024b and earlier
    end
end
