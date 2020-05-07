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
                         'interpFace', [], ...
                         'OSflux', []);
            opt = merge_options(opt, varargin{:});

            % Set up HAP
            if isempty(opt.interpFace)
                ntpfa.interpFace = findHAP(model.G, model.rock);
                ntpfa.interpFace = correctHAP(model.G, ntpfa.interpFace, opt.myRatio);
            else
                ntpfa.interpFace = opt.interpFace;
            end
            dispif(mrstVerbose, 'fraction of faces with HAPs outside convex hull is %d\n', ntpfa.interpFace.fraction);

            % Set up one-sided fluxes
            if isempty(opt.OSflux)
                ntpfa.OSflux = findOSflux(model.G, model.rock, ntpfa.interpFace);
            else
                ntpfa.OSflux = opt.OSflux;
            end

            ntpfa.scale = -1 ./ model.operators.T;
            ntpfa.internalConn = model.operators.internalConn;
            ntpfa.nph = sum(model.getActivePhases());
        end

        function v = gradient(ntpfa, pressure)
            grad = cell(1, ntpfa.nph);
            for i = 1:ntpfa.nph
                T = TransNTPFA(ntpfa.G, pressure, ntpfa.OSflux);
                v = computeFlux(ntpfa.G, pressure, T);
                v = v(ntpfa.internalConn, :);
                v = v .* ntpfa.scale;
                grad{i} = v;
            end
        end
    end

end
