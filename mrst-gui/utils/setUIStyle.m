function [] = setUIStyle(d, styleName, layoutOpt)
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
                    'BorderWidth',            0);
                
dummyPanelProps = panelProps;                
   
% set default parameters
if isa(d, 'uiItem')
    d.layout.params = getDefaultsItemParams();
elseif isa(d, 'uiMenu')
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
        if isa(d, 'uiMenu') && d.level == 1
            d.layout.params.outerMargins    = [1 1 1 1];
        end
        if isa(d, 'uiMenu')
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
        if isa(d, 'uiItem')
            props.titleBackgroundColor = blue;
            props.BackgroundColor = standard;
            props.titleColor = standard;
        else
            props.BackgroundColor = blue;
            props.ForegroundColor = standard;
        end
        if isa(d, 'uiMenu') && d.level == 1
            d.layout.params.dummyOuterMargins    = [2 2 2 2];
        else
            d.layout.params.dummyOuterMargins    = [2 2 2 2];
        end
        if isa(d, 'uiMenu')
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
        if isa(d, 'uiMenu')
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
if isa(d, 'uiMenu')
    d.titleHeight = ceil(18*d.FontSize/10) + yMrg + panelProps.BorderWidth;
else
    d.titleHeight = ceil(18*d.titleFontSize/10) + yMrg + panelProps.BorderWidth;
end   
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
    params = struct('margins',            [5 0 0 0], ...  % left, right, bottom, top
                    'outerMargins',       [0 0 0 0], ...  % left, right, bottom, top 
                    'dummyOuterMargins',  [0 0 0 0], ...  %
                    'vskipMenu',                  1, ...
                    'vskipItem',                  5);       
end

