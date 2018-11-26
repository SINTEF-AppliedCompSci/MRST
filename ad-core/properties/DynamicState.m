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

        function state = getState(dyn_state)
            state = dyn_state.state;
        end
        
        function u = subsasgn(u, s, v)
            % Subscripted reference. Called for `u(s) = v`
            if strcmp(s.type, '.') && ischar(s.subs)
                act = strcmp(s.subs, u.dynamicNames);
                if any(act)
                    u.dynamicVariables{act} = v;
                else
                    u.state.(s.subs) = v;
                end
            else
                u = builtin('subsasgn', u, s,v);
            end
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