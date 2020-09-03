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
