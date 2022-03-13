function setup = spatialUpscalingConstructor(setup, partition, varargin)
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

    opt = struct('cartBox', false);
    opt = merge_options(opt, varargin{:});

    nc = setup.model.G.cells.num;

    model = upscaleModelTPFA(setup.model, partition, 'validatePartition', false);
    if opt.cartBox
        model.G.cartDims = ones(1, model.G.physDims)*sqrt(model.G.cells.num, 1);
    end
    % Map well cells
    schedule = setup.schedule;
    for i = 1:numel(schedule.control)
        W = setup.schedule.control(i).W;
        if ~isempty(W)
            for j = 1:numel(W)
               W(j).cells = model.G.partition(W(j).cells);
            end
            schedule.control(i).W = W;
        end
    end
    % Map state
    fnames = fieldnames(setup.state0);
    donor  = setup.state0;
    state0 = struct();
    for i = 1:numel(fnames)
        n = fnames{i};
        if isa(donor.(n), 'double') && size(donor.(n), 1) == nc
            state0.(n) = zeros(model.G.cells.num, size(donor.(n),2));
            state0.(n)(model.G.partition,:) = donor.(n);
        end
    end

    setup.state0 = state0;
    setup.model = model;
    setup.schedule = schedule;

end
