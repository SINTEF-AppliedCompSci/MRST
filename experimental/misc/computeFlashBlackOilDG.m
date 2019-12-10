function state = computeFlashBlackOilDG(state, state0, model, status)

    statePrev = state;
    % Flash
    state = computeFlashBlackOil(state, state0, model, status);
    names = {'s', 'rs', 'rv'};
    for vNo = 1:numel(names)
        v = names{vNo};
        if isfield(state, v) && size(state.(v),1) == model.G.cells.num
            state = scaleDof(model.disc, state, statePrev, v);
        end
    end

end

function state = scaleDof(disc, state, statePrev, v)

    vdof = [v, 'dof'];
    % Add to change to constant part
    d  = state.(v) - statePrev.(v);
    ix = disc.getDofIx(state, 1, Inf);
    state.(vdof)(ix,:) = state.(vdof)(ix,:) + d;
    % Adjust mean
    for cNo = 1:size(state.(vdof),2)
        val = disc.getCellMean(state, state.(vdof)(:,cNo));
        f   = state.(v)(:,cNo)./val;
        f(~isfinite(f)) = 1;
        state.(vdof)(:,cNo) = state.(vdof)(:,cNo).*rldecode(f, state.nDof, 1);
    end
    
    if 1
        for cNo = 1:size(state.(vdof),2)
            val = disc.getCellMean(state, state.(vdof)(:,cNo));
            d   = val - state.(v)(:,cNo);
        end
    end

end