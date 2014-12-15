function [sol, withinLims] = updateSwitchedWellControls(wellmodel, model, sol, pBH, q_s)
W = wellmodel.W;
if isempty(W)||isempty(W(1).lims)
    withinLims = true(numel(sol),1);
else
    allowWellSignChange = wellmodel.allowWellSignChange;

    nwells     = numel(sol);
    withinLims = true(nwells,1);

    pBH   = double(pBH);
    q_s   = cell2mat( cellfun(@double, q_s, 'UniformOutput', false) );

    for wnr = 1:numel(sol)
        if isfield(W(wnr), 'status') && ~W(wnr).status
            % Inactive well, skip any limit checks
            continue
        end
        lims = W(wnr).lims;
        if ~allowWellSignChange
            lims.vrat = -inf;
        else
            lims.vrat = -inf;
        end
        pBHw  = pBH(wnr);
        q_sw  = q_s(wnr,:);
        qt_sw = sum(q_sw);
        if ~isnumeric(W(wnr).lims)
            if sol(wnr).sign > 0   % injector
                modes   = {'bhp', 'rate', 'rate'};
                flags = [pBHw > lims.bhp, qt_sw > lims.rate, qt_sw < lims.vrat];
            else            % producer
                modes   = {'bhp', 'orat', 'lrat', 'grat', 'wrat', 'vrat'};
                flags = [pBHw       < lims.bhp,  ...
                    q_sw(2)         < lims.orat, ...
                    q_sw(1)+q_sw(2) < lims.lrat, ...
                    q_sw(3)         < lims.grat, ...
                    q_sw(1)         < lims.wrat, ...
                    qt_sw           > -lims.vrat];
            end
        else
            modes = {};
            flags = false(numel(sol), 1);
            assert(isinf(lims))
        end
        %limits we need to check (all others than w.type):
        chkInx = ~strcmp(sol(wnr).type, modes);
        vltInx = find(flags(chkInx), 1);
        if ~isempty(vltInx)
            withinLims(wnr) = false;
            modes  = modes(chkInx);
            switchMode = modes{vltInx};
            fprintf('Well %s: Control mode changed from %s to %s.\n', sol(wnr).name, sol(wnr).type, switchMode);
            sol(wnr).type = switchMode;
            sol(wnr).val  = lims.(switchMode);
        end
    end
end
end


