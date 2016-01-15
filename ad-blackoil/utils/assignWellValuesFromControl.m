function wellSol = assignWellValuesFromControl(model, wellSol, W, wi, oi, gi)
% Utility function to assign well values that are the resut of the controls
% being enabled.
    for w = 1:numel(wellSol)
        ws = wellSol(w);
        tp = ws.type;
        v  = ws.val;
        switch tp
            case 'bhp'
                ws.bhp = v;
            case 'rate'
                if model.water
                    ws.qWs = v*W(w).compi(wi);
                end
                if model.oil
                    ws.qOs = v*W(w).compi(oi);
                end
                if model.gas
                    ws.qGs = v*W(w).compi(gi);
                end
                if isprop(model, 'polymer') && model.polymer
                    ws.qWPoly = ws.qWs*W(w).poly;
                end
                if isprop(model, 'surfactant') && model.surfactant
                    ws.qWSurfact = ws.qWs*W(w).surfact;
                end
            case 'orat'
                ws.qOs = v;
            case 'wrat'
                ws.qWs = v;
            case 'grat'
                ws.qGs = v;
            case 'lrat'
                % Do nothing
            otherwise
                error('Unknown well control mode');
        end
        wellSol(w) = ws;
    end
end