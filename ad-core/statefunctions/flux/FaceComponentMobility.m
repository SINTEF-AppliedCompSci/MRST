classdef FaceComponentMobility < StateFunction & UpwindProperty
    properties (Access = protected)
        upwind_name; % Name of state function where upwind flag comes from
    end
    
    methods
        function fm = FaceComponentMobility(model, upwinding, upwind_name)
            if nargin < 3
                upwind_name = 'PhaseUpwindFlag';
            end
            fm@StateFunction(model);
            fm@UpwindProperty(upwinding)
            fm.upwind_name = upwind_name;
            fm = fm.dependsOn(upwind_name);
            fm = fm.dependsOn('ComponentMobility', 'FlowPropertyFunctions');
        end
        
        function mobf = evaluateOnDomain(prop, model, state)
            flag = prop.getEvaluatedDependencies(state, prop.upwind_name);
            mob = model.getProps(state, 'ComponentMobility');
            [ncomp, nph] = size(mob);
            mobf = cell(ncomp, nph);
            for c = 1:ncomp
                for ph = 1:nph
                    if ~isempty(mob{c, ph})
                        mobf{c, ph} = prop.faceUpstream(model, state, flag{ph}, mob{c, ph});
                    end
                end
            end
        end
    end
end

%{
Copyright 2009-2019 SINTEF Digital, Mathematics & Cybernetics.

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
