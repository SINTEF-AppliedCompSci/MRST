classdef ComponentProperty
    % Virtual class for properties that get their values from the model's
    % components
    properties (Access = protected)
        componentFunctionName % Name of member function in component implementation
        includeStandard = true;
        includeShrinkageFactors = [];
    end
    
    methods
        function gp = ComponentProperty(model, name)
            if nargin > 0
                gp.componentFunctionName = name;
                ncomp = model.getNumberOfComponents();
                deps = cell(ncomp, 1);
                exts = cell(ncomp, 1);
                cfn = gp.componentFunctionName;
                for c = 1:ncomp
                    cdeps = model.Components{c}.getFunctionDependencies(cfn);
                    deps{c} = cdeps.dependencies;
                    exts{c} = cdeps.externals;
                end
                % Internal - flow props dependencies
                deps = unique(vertcat(deps{:}));
                gp = gp.dependsOn(deps); %#ok virtual class
                % Manage external dependencies
                exts = vertcat(exts{~cellfun(@isempty, exts)});
                if ~isempty(exts)
                    names = {exts.name};
                    [~, pos] = unique(names);
                    exts = exts(pos);
                    gp = gp.dependsOn(exts); %#ok virtual class
                end
                if isempty(gp.includeShrinkageFactors)
                    gp.includeShrinkageFactors = isa(model, 'ThreePhaseBlackOilModel');
                end
                if gp.includeStandard
                    if gp.includeShrinkageFactors
                        gp = gp.dependsOn({'ShrinkageFactors'}, 'PVTPropertyFunctions');
                    end
                    gp = gp.dependsOn({'Density', 'PoreVolume'}, 'PVTPropertyFunctions');
                    gp = gp.dependsOn('Mobility', 'FlowPropertyFunctions');
                end
            else
                gp.includeShrinkageFactors = false;
            end
        end

        function extra = getExtraArguments(prop, model, state)
            if prop.includeStandard
                [rho, pv, mob] = prop.getEvaluatedExternals(model, state, 'Density', 'PoreVolume', 'Mobility');
                s = model.getProp(state, 's');
                extra = struct('s', [], 'pv', pv, 'rho', [], 'mob', []);
                extra.s = s;
                extra.rho = rho;
                extra.mob = mob;
                if prop.includeShrinkageFactors
                    extra.b = prop.getEvaluatedExternals(model, state, 'ShrinkageFactors');
                end
                extra = {extra};
            else
                extra = {};
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
