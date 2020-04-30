classdef TwoPhaseNPVSelector < UIItem
    properties
        Callback
        roEdit
        rwpEdit
        rwiEdit
        discountEdit
    end
    properties (Dependent)
        ro
        rwp
        rwi
        d
    end
    
    methods
        
        function s = TwoPhaseNPVSelector(varargin)
            opt = struct('Parent',          [], ...
                         'Callback',        @(src, event)disp('Boo'), ...
                         'Position',        [10 10 300 200], ...
                         'Visible',         'on', ...
                         'Title',           'NPV options', ...
                         'ro',              60, ...
                         'rwp',              5, ...
                         'rwi',             10, ...
                         'discount',         0);
            [opt, extraOpt] = merge_options(opt, varargin{:});
            
            textOpts = {'Parent', [], 'Style', 'text',    'Value',  [], 'Visible', 'off'};
            editOpts = {'Parent', [], 'Style', 'edit',    'Value',  [], 'Visible', 'off'};
            
            roText   = uicontrol(textOpts{:}, 'String', 'Oil revenue [$/bbl]:');
            roEdit   = uicontrol(editOpts{:}, 'String', num2str(opt.ro));
            
            rwpText  = uicontrol(textOpts{:}, 'String', 'Water production cost [$/bbl]');
            rwpEdit  = uicontrol(editOpts{:}, 'String', num2str(opt.rwp));
            
            rwiText  = uicontrol(textOpts{:}, 'String', 'Water injection cost [$/bbl]');
            rwiEdit  = uicontrol(editOpts{:}, 'String', num2str(opt.rwi));
            
            discountText  = uicontrol(textOpts{:}, 'String', 'Discount factor [%/year]');
            discountEdit  = uicontrol(editOpts{:}, 'String', num2str(opt.discount));     

                               
            controls      = {{roText,  roEdit}, ...
                             {rwpText, rwpEdit}, ...
                             {rwiText, rwiEdit}, ...
                             {discountText, discountEdit}};

            controlLayout = {[.75 .25], ...
                             [.75 .25], ...
                             [.75 .25], ...
                             [.75 .25]}; 
      
            s = s@UIItem('Parent', opt.Parent, 'controls', controls, ...
                                'controlWidths', controlLayout, 'Title', opt.Title, ...
                                'Position', opt.Position,'Visible', 'off', extraOpt{:});
            
            s.roEdit        = roEdit;
            s.rwpEdit       = rwpEdit;
            s.rwiEdit       = rwiEdit;
            s.discountEdit  = discountEdit;   
            s.fixedHeight   = true;
            
            % main callback
            s.Callback = opt.Callback;
            % all edits have same callback
            set([s.roEdit, s.rwpEdit, s.rwiEdit, s.discountEdit], 'Callback', @s.editCallback);
            set([s.roEdit, s.rwpEdit, s.rwiEdit, s.discountEdit], 'HorizontalAlignment', 'right');
            % set visible
            s.Visible = opt.Visible;
        end

        function set.ro(s, val)
            val = capValue(val, [0 1000]);
            if isfinite(val)
                s.roEdit.String = num2str(val);
            end
        end
        function val = get.ro(s)
            val = str2double(s.roEdit.String);
        end
        
        function set.rwp(s, val)
            val = capValue(val, [0 1000]);
            if isfinite(val)
                s.rwpEdit.String = num2str(val);
            end
        end
        function val = get.rwp(s)
            val = str2double(s.rwpEdit.String);
        end
        
        function set.rwi(s, val)
            val = capValue(val, [0, 1000]);
            if isfinite(val)
                s.rwiEdit.String = num2str(val);
            end
        end
        function val = get.rwi(s)
            val = str2double(s.rwiEdit.String);
        end
        
        function set.d(s, val)
            val = capValue(val, [0, 99]);
            if isfinite(val)
                s.discountEdit.String = num2str(val);
            end
        end
        function val = get.d(s)
            val = str2double(s.discountEdit.String);
        end
        
        function editCallback(s, src, event)
            s.Callback(src, event)
        end
    end
end

function v = capValue(v, lims)
if ischar(v)
    v = str2double(v);
end
assert(isnumeric(v), 'Non-numeric value...')
v = min(lims(2), max(lims(1), v));
end
