classdef DynamicState 
    properties
        state
        dynamicNames = {};
        dynamicVariables = {};
        implicitEvaluation = true
    end
    
    methods
        function dyn = DynamicState(state, varargin)
            % Calling syntax:
            % DynamicState(state)
            % DynamicState(state, names, variables)
            % DynamicState(..., varargin)
            if nargin > 1
                if ~ischar(varargin{1})
                    dyn.dynamicNames = varargin{1};
                    dyn.dynamicVariables = varargin{2};
                end
            end
            dyn.state = state;
            dyn.implicitEvaluation = numel(dyn.dynamicVariables) > 0;
        end
        
        function addDynamicVariables(variables, names)
            
        end
        
        function state = getState(dyn_state)
            state = dyn_state.state;
        end
        
        function h = subsref(u, s)
            % h = u(s)
            fld = s(1).subs;
            if strcmp(s(1).type, '.') && ischar(fld)
                if strcmp(s(1).type, 'state')
                    h = u.state;
                    if numel(s) > 1
                        h = builtin('subsref', h, s(2:end));
                    end
                    return
                end
                
                act = strcmp(u.dynamicNames, fld);
                if any(act)
                    dv = u.dynamicVariables{act};
                    if numel(s) > 1
                        if iscell(dv) && strcmp(s(2).type, '()')
                            s(2).type = '{}';
                        end
                        h = builtin('subsref', dv, s(2:end));
                    else
                        h = dv;
                    end
                elseif isfield(u.state, fld)
                    h = builtin('subsref', u.state, s);
                else
                    h = builtin('subsref', u, s);
                end
            else
                h = builtin('subsref',u,s);
            end
        end
    end
end