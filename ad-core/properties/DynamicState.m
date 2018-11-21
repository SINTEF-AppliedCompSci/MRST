classdef DynamicState < handle
    properties
        state
        dynamicNames = {};
        dynamicVariables = {};
        evaluatedProperties = struct();
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
        
        function h = subsref(u, s)
            % h = u(s)
            if strcmp(s(1).type, '.') && ischar(s(1).subs)
                act = strcmp(u.dynamicNames, s(1).subs);
                if any(act)
                    dv = u.dynamicVariables{act};
                    if numel(s) > 1
                        h = builtin('subsref', dv, s(2:end));
                    else
                        h = dv;
                    end
                else
                    h = builtin('subsref', u.state, s);
                end
            else
                h = builtin('subsref',u,s);
            end
        end
    end
end