classdef Hybrid
    properties

        G
        fdmodels
        faceblocks
        op
        facemap

    end

    methods

        function hybrid = Hybrid(model, fdmodels, faceblocks, varargin)

            hybrid.G = model.G;
            hybrid.fdmodels = fdmodels;
            hybrid.faceblocks = faceblocks;
            hybrid.op = model.operators;

            G = hybrid.G;
            op = hybrid.op;
            internal = find(op.internalConn);
            hybrid.facemap = sparse(internal, 1:numel(internal), 1, G.faces.num, numel(internal));

        end

        function grad = gradient(hybrid, pressure)

            fdmodels = hybrid.fdmodels;
            faceblocks = hybrid.faceblocks;
            op = hybrid.op;
            facemap = hybrid.facemap;

            g = fdmodels{1}.PressureGradient.Grad(pressure);
            grad = facemap*g;

            for k = 2:numel(fdmodels)
                g = fdmodels{k}.PressureGradient.Grad(pressure);
                g = facemap*g;
                grad(faceblocks{k}) = g(faceblocks{k});
            end

            grad = grad(op.internalConn);

        end

    end

end

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
