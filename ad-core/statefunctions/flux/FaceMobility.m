classdef FaceMobility < StateFunction & UpwindProperty
    % Phase mobility, with phase upwind (or alternatively some other
    % upwind, provided that the name of that flag is given to the
    % constructor)
    properties (Access = protected)
        upwind_name; % Name of state function where upwind flag comes from
    end
    
    methods
        function fm = FaceMobility(model, upstr, upwind_name)
            if nargin < 2
                upstr = UpwindFunctionWrapperDiscretization(model);
            end
            if nargin < 3
                upwind_name = 'PhaseUpwindFlag';
            end
            fm@StateFunction(model);
            fm@UpwindProperty(upstr);
            fm.upwind_name = upwind_name;
            fm = fm.dependsOn(upwind_name);
            fm = fm.dependsOn('Mobility', 'FlowPropertyFunctions');
            fm.label = '\lambda_\alpha^f';
        end
        
        function fmob = evaluateOnDomain(prop, model, state)
            flag = prop.getEvaluatedDependencies(state, prop.upwind_name);
            mob = prop.getEvaluatedExternals(model, state, 'Mobility');
            nph = numel(mob);
            fmob = cell(1, nph);
            for i = 1:nph
                fmob{i} = prop.faceUpstream(model, state, flag{i}, mob{i});
            end
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
