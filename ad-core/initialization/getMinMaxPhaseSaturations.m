function [s_min, s_max] = getMinMaxPhaseSaturations(model, cellInx)
    if nargin < 2
        cellInx = 1;
    end
    nc = numel(cellInx);
    [s_min, s_max] = getMinMaxPhaseSaturationsFromRelPerm(model, 1e-6, cellInx(1));
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
                s_min(:, wix) = f.sWcon(cellInx);
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
