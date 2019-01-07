function [s_min, s_max] = getMinMaxPhaseSaturations(model, satnum)
    if nargin < 2
        satnum = 1;
    end
    nc = numel(satnum);
    if isfield(model.fluid, 'krPts')
        pts = model.fluid.krPts;
        phases = model.getPhaseNames();
        nph = numel(phases);
        [s_min, s_max] = deal(1, nph);
        for i = 1:nph
            ph = lower(phases(i));
            if strcmp(ph, 'o') && ~isfield(pts, 'o')
                if model.water && model.gas
                    s_min(i) = min(pts.ow(satnum, 2), pts.og(satnum, 2));
                    s_max(i) = max(pts.ow(satnum, 3), pts.og(satnum, 3));
                elseif model.water
                    s_min(i) = pts.ow(satnum, 2);
                    s_max(i) = pts.ow(satnum, 3);
                elseif model.gas
                    s_min(i) = pts.og(satnum, 2);
                    s_max(i) = pts.og(satnum, 3);
                end
            else
                s_min(i) = pts.(ph)(satnum, 2);
                s_max(i) = pts.(ph)(satnum, 3);
            end
        end
    else
        [s_min, s_max] = getMinMaxPhaseSaturationsFromRelPerm(model, 1e-6, cellInx(1));
    end

    s_min = repmat(s_min, nc, 1);
    s_max = repmat(s_max, nc, 1);
    f = model.fluid;
    wix = model.getPhaseIndex('W');
    oix = model.getPhaseIndex('O');
    gix = model.getPhaseIndex('G');
    if model.gas
        % Minimum gas should always be zero
        s_min(:, gix) = 0;
    end
    
    if model.water
        if isfield(f, 'sWcon')
            % Minimum water = connate water
            if numel(f.sWcon) == 1
                s_min(:, wix) = f.sWcon;
            else
                s_min(:, wix) = f.sWcon(satnum);
            end
            % Water can always be one
            s_max(:, wix) = 1;
        end
        if model.gas
            % Account for minimum water saturation
            s_max(:, gix) = min(s_max(:, gix), 1 - s_min(:, wix));
        end
    end
end
