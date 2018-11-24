function [s_min, s_max] = getMinMaxPhaseSaturations(model, satnum)
    if nargin < 2
        satnum = 1;
    end
    nc = numel(satnum);
    if isfield(model.fluid, 'krPts')
        pts = model.fluid.krPts;
        s_min = [pts.w(satnum, 2), max(pts.ow(satnum, 2), pts.og(satnum, 2)), pts.g(satnum, 2)];
        s_max = [pts.w(satnum, 3), max(pts.ow(satnum, 3), pts.og(satnum, 3)), pts.g(satnum, 3)];
    else
        [s_min, s_max] = getMinMaxPhaseSaturationsFromRelPerm(model, 1e-6, cellInx(1));
    end
%     s_max = ones(size(s_min));

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
        end
        if model.gas
            % Account for minimum water saturation
            s_max(:, gix) = min(s_max(:, gix), 1 - s_min(:, wix));
        end
    end
end
