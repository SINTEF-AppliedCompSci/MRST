classdef NTPFA
    properties
        interpFace % Harmonic averaging points
        OSflux % One side fluxes
        bc % homogeneous Neumann BCs (Kobaisi format)
        G
        scale
        internalConn
        nph
    end

    methods

        function ntpfa = NTPFA(model, opt)
            ntpfa.G = model.G;

            % Set hom Neumann bc
            ntpfa.bc.face = boundaryFaces(model.G);
            ntpfa.bc.type = repmat({'flux'}, numel(ntpfa.bc.face), 1);
            ntpfa.bc.value = repmat({@(x) 0}, numel(ntpfa.bc.face), 1);

            % Set up HAP and fluxes
            ntpfa.interpFace = findHAP(model.G, model.rock, ntpfa.bc);
            ntpfa.interpFace = correctHAP(model.G, ntpfa.interpFace, opt.myRatio);
            disp(['fraction of faces with HAPs outside convex hull ', ...
                num2str(ntpfa.interpFace.fraction)])
            ntpfa.OSflux = findOSflux(model.G, model.rock, ntpfa.bc, ntpfa.interpFace);

            ntpfa.scale = -1 ./ model.operators.T;
            ntpfa.internalConn = model.operators.internalConn;
            ntpfa.nph = sum(model.getActivePhases());
        end

        function v = gradient(ntpfa, pressure)
            grad = cell(1, ntpfa.nph);
            for i = 1:ntpfa.nph
                T = TransNTPFA(ntpfa.G, ntpfa.bc, [], pressure, ntpfa.OSflux);
                v = computeFlux(ntpfa.G, pressure, T, [], [], []);
                v = v(ntpfa.internalConn, :);
                v = v .* ntpfa.scale;
                grad{i} = v;
            end
        end
    end
    
end