classdef AvgMPFA < NTPFA
% Average AvgMPFA

    properties

        % T is a cell structure containing two matrices which compute the fluxes. Each
        % matrix maps cell pressure values to internal face flux value. We have
        % one matrix for each side of the face (the first matrix for
        % G.faces.neighbors(: , 1) and the second for G.faces.neighbors(:, 2)).

        T

    end

    methods

        function avgmpfa = AvgMPFA(model, varargin)

            avgmpfa = avgmpfa@NTPFA(model, varargin{:});
            G = model.G;
            OSflux = avgmpfa.OSflux;
            avgmpfa.T = AvgTransNTPFA(G, OSflux);

        end

        function grad = gradient(ntpfa, pressure)

            grad = cell(1, ntpfa.nph);
            T = ntpfa.T;

            grad = 0.5*(T{1}*pressure - T{2}*pressure);
            grad = grad .* ntpfa.scale;

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
