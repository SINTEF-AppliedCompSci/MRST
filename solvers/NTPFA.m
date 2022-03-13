classdef NTPFA
    properties
        interpFace % Harmonic averaging points
        OSflux % One-side fluxes
        G
        scale
        internalConn
        nph
        lagged
    end

    methods

        function ntpfa = NTPFA(model, varargin)
            ntpfa.G = model.G;

            opt = struct('myRatio', [], ...
                         'OSflux', [], ...
                         'interpFace', [], ...
                         'saveOSfluxDir', [], ...
                         'saveInterpFaceDir', [], ...
                         'lagged', false);

            opt = merge_options(opt, varargin{:});

            ntpfa.lagged = opt.lagged;

            % Set up HAP
            if isempty(opt.interpFace)
                interpFace = findHAP(model.G, model.rock);
                dispif(mrstVerbose, ...
                       'fraction of cells with centroids outside convex hull is %d\n', ...
                       interpFace.fraction);

                interpFace = correctHAP(model.G, interpFace, opt.myRatio);

                if ~isempty(opt.saveInterpFaceDir)
                    save(fullfile(opt.saveInterpFaceDir, 'interpFace.mat'), ...
                         'interpFace');
                end

                ntpfa.interpFace = interpFace;
            else
                ntpfa.interpFace = opt.interpFace;
                dispif(mrstVerbose, ...
                       'fraction of faces with HAPs outside convex hull is %d\n', ...
                       ntpfa.interpFace.fraction);
            end

            % Set up one-sided fluxes
            if isempty(opt.OSflux)
                OSflux = findOSflux(model.G, model.rock, ntpfa.interpFace);

                if ~isempty(opt.saveOSfluxDir)
                    save(fullfile(opt.saveOSfluxDir, 'OSflux.mat'), ...
                         'OSflux');
                end

                ntpfa.OSflux = OSflux;
            else
                ntpfa.OSflux = opt.OSflux;
            end

            ntpfa.scale = 1./model.operators.T;
            ntpfa.internalConn = model.operators.internalConn;
            if ismethod(model, 'getActivePhases')
                ntpfa.nph = sum(model.getActivePhases());
            end
        end

        function grad = gradient(ntpfa, pressure)
            grad = ntpfa.unscaledGradient(pressure);
            grad = grad .* ntpfa.scale;
        end

        function grad = unscaledGradient(ntpfa, pressure)

            if ntpfa.lagged
                T = TransNTPFA(ntpfa.G, value(pressure), ntpfa.OSflux);
            else
                T = TransNTPFA(ntpfa.G, pressure, ntpfa.OSflux);
            end

            v = computeFlux(ntpfa.G, pressure, T);
            v = v(ntpfa.internalConn, :);
            grad = -v;
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
