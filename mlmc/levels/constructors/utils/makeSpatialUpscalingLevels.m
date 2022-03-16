function levels = makeSpatialUpscalingLevels(setup, varargin)
%Undocumented Utility Function

%{
Copyright 2009-2022 SINTEF Digital, Mathematics & Cybernetics.

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

    opt = struct('numLevels', 3);
    [opt, extra] = merge_options(opt, varargin{:});
    model = setup.model;
    partitions0 = makeNestedPartitions2(model, 2*opt.numLevels-1, extra{:});
    partitions0 = partitions0(1:2:end);
    partitions  = partitions0;
    for i = 2:opt.numLevels
        level = repmat(i, setup.model.G.cells.num, 1);
        level(vertcat(setup.schedule.control(1).W.cells)) = 1;
        partitions{i} = makeHybridPartition(partitions0, level, inf, []);
    end
    partitions = flipud(partitions);
    
    levels = cell(opt.numLevels, 1);
    for i = 1:opt.numLevels
        levels{i} = SpatialUpscalingLevel(i, partitions{i});
    end
end
