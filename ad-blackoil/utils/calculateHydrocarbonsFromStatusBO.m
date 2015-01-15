function [sG, rs, rv, rsSat, rvSat] = calculateHydrocarbonsFromStatusBO(model, status, sO, x, rs, rv, pressure)
    % define sG, rs and rv in terms of x
    % We cast "sO" to double to avoid cancellation of derivatives when
    % calculating the final sO later.
    sG = status{2}.*double(sO) + status{3}.*x;
    if model.disgas
        rsSat = model.fluid.rsSat(pressure);
        rs = (~status{1}).*rsSat + status{1}.*x;
    else % otherwise rs = rsSat = const
        rsSat = rs;
    end
    if model.vapoil
        rvSat = model.fluid.rvSat(pressure);
        rv = (~status{2}).*rvSat + status{2}.*x;
    else % otherwise rv = rvSat = const
        rvSat = rv;
    end
end
