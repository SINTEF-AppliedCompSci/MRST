classdef PreDiagnosticsSelector < UIItem
    properties
        selector
        intervalPopup
        minSlider
        minText
        maxSlider
        maxText
        OKButton
        t
        t_date
    end
    properties (Dependent, SetObservable)
        Value = []
        ix
    end
    
    methods
        
        function s = preDiagnosticsSelector(varargin)
            
            opt = struct('Parent',          [], ...
                         'Visible',         'on', ...
                         'Title', 'Select time-steps for diagnostics', ...
                         'restartInfo', []);
            [opt, extraOpt] = merge_options(opt, varargin{:});
                        
            % get times in days:
            assert(~isempty(opt.restartInfo), 'Need restart info to select time steps ...')
            info = opt.restartInfo;
            startday = datenum(info.date(1, [3 2 1]));
            t        = startday + info.time - info.time(1);
            n = numel(t);
            t_date   = cellfun(@(x)datestr(x, 'mmm dd, yyyy'), mat2cell(t, ones(1, n)), 'UniformOutput', false); 
            
            selector  = uicontrol('Parent', [], 'Style', 'listbox',  'Max', 2, 'Min', 0, ...
                                       'Value', [], 'String', t_date, 'Visible', 'off');
            minSlider = uicontrol('Parent', [], 'Style', 'slider',...
                                  'Min', t(1), 'Max', t(end), 'Value', t(1));
            minText   = uicontrol('Parent', [], 'Style', 'text', 'Value', [], ... 
                                      'String', t_date{1}, 'Visible', 'off');
            maxSlider = uicontrol('Parent', [], 'Style', 'slider',...
                                  'Min', t(1), 'Max', t(end), 'Value', t(end));
            maxText   = uicontrol('Parent', [], 'Style', 'text', 'Value', [], ... 
                                      'String', t_date{end}, 'Visible', 'off');
            
            intText   = uicontrol('Parent', [], 'Style', 'text', 'Value', [], ... 
                                      'String', 'Select time-steps at (approximate) intervals:', 'Visible', 'off');
            intervalPopup     = uicontrol('Parent', [], 'Style', 'popup', ...
                                          'String', {'n/a', 'month', '6 months', 'year', '5 year'},...
                                          'Value', 1, 'UserData', [nan 30 365/2 365 5*365]);
                                      
            OKButton          = uicontrol('Parent', [], 'Style', 'pushbutton', 'Tag', 'play', 'String', 'OK');    
             
            controls = {{selector}, {minSlider, minText}, {maxSlider, maxText}, {intText}, {intervalPopup}, {[], OKButton}};
            layout   = {nan, [.7 .3], [.7 .3], nan, nan, [nan 0.3]};
            pos = [0 0 opt.Parent.Position(3:4)];
            s = s@uiItem('Parent', opt.Parent, 'controls', controls, 'controlWidths', layout, ...
                               'Title', opt.Title, 'Position', pos,'Visible', 'on', extraOpt{:});
                                   
            s.selector      = selector;
            s.intervalPopup = intervalPopup;
            s.minSlider     = minSlider;
            s.minText       = minText;
            s.maxSlider     = maxSlider;
            s.maxText       = maxText;
            s.OKButton      = OKButton;
            s.t = t;
            s.t_date = t_date;
            % main callback
            selector.Callback      = @s.selectorCallback;
            intervalPopup.Callback = @s.intervalCallback;
            minSlider.Callback     = @s.minCallback;
            maxSlider.Callback     = @s.maxCallback;
            OKButton.Callback      = @s.buttonCallback;
            % set visible
            s.Visible = opt.Visible;
            
            % add item context menues after figure has been created
            s.selector.UIContextMenu  = s.listboxContextMenu({'Select all steps', 'Clear selected steps'}, ...
                                                               {@s.selectAll, @s.clearAll} ); 
        end
        
        function set.Value(s, val)
            nv = numel(s.selector.String);
            assert(all(val>0) && all(val<=nv));
            s.selector.Value = val;
        end
        function val = get.Value(s)
            val = s.selector.Value;
        end
        
        function set.ix(s, val)
            s.Value = val;
        end
        function val = get.ix(s)
            val = s.Value;
        end
        
        function selectAll(s, src, event)
            nv = numel(s.selector.String);
            if nv ~= numel(s.Value)
                s.Value = (1:nv);
                s.selectorCallback(src, event)
            end
        end
        
        function clearAll(s, src, event)
            if ~isempty(s.Value)
                s.Value = [];
                s.selectorCallback(src, event)
            end
        end
        
        function selectorCallback(s, src, event)
            s.ix = s.selector.Value;
            if ~isempty(s.ix)
                range = s.t( s.ix([1 end]) );
                if s.minSlider.Value > range(1)
                     s.minCallback(src, event);
                elseif s.maxSlider.Value  < range(2)
                    s.maxCallback(src, event);
                end
            end
        end
        
        function intervalCallback(s, src, event)
            if ~(s.intervalPopup.Value == 1)
                range = [s.minSlider.Value, s.maxSlider.Value];
                dt = s.intervalPopup.UserData(s.intervalPopup.Value);
                sel = range(1):dt:range(2);
                % find corresponing indices
                ii = zeros(size(sel));
                for k = 1:numel(ii)
                    [~, ii(k)] = min(abs(s.t-sel(k)));
                end
                s.ix = unique(ii);
            end
        end
        
        function minCallback(s, src, event)
            s.minSlider.Value = min(s.minSlider.Value, s.maxSlider.Value);
            keep = s.t(s.ix) >= s.minSlider.Value;
            s.ix = s.ix(keep);
            s.minText.String = datestr(s.minSlider.Value, 'mmm dd, yyyy');
            s.intervalCallback(src, event);
        end
        
        function maxCallback(s, src, event)
            s.maxSlider.Value = max(s.minSlider.Value, s.maxSlider.Value);
            keep = s.t(s.ix) <= s.maxSlider.Value;
            s.ix = s.ix(keep);
            s.maxText.String = datestr(s.maxSlider.Value, 'mmm dd, yyyy');
            s.intervalCallback(src, event);
        end
        
        function buttonCallback(s, src, event)
            % call close function for parent
            feval(s.Parent.CloseRequestFcn);
        end
    end
end
