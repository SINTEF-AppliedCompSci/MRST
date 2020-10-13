classdef NTPFA
    properties
        interpFace % Harmonic averaging points
        OSflux % One-side fluxes
        G
        scale
        internalConn
        nph
    end

    methods

        function ntpfa = NTPFA(model, varargin)
            ntpfa.G = model.G;

            opt = struct('myRatio', [], ...
                         'OSflux', [], ...
                         'interpFace', [], ...
                         'saveOSfluxDir', [], ...
                         'saveInterpFaceDir', []);

            opt = merge_options(opt, varargin{:});

            % Set up HAP
            if isempty(opt.interpFace)
                interpFace = findHAP(model.G, model.rock);
                interpFace = correctHAP(model.G, interpFace, opt.myRatio);

                if ~isempty(opt.saveInterpFaceDir)
                    save(fullfile(opt.saveInterpFaceDir, 'interpFace.mat'), ...
                         'interpFace');
                end

                ntpfa.interpFace = interpFace;
            else
                ntpfa.interpFace = opt.interpFace;
            end
            dispif(mrstVerbose, ...
                   'fraction of faces with HAPs outside convex hull is %d\n', ...
                   ntpfa.interpFace.fraction);

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
            T = TransNTPFA(ntpfa.G, pressure, ntpfa.OSflux);
            v = computeFlux(ntpfa.G, pressure, T);
            v = v(ntpfa.internalConn, :);
            grad = -v;
        end
    end

end

%{
Copyright 2009-2020 SINTEF Digital, Mathematics & Cybernetics.

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
