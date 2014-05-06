function [sol, withinLims] = updateControls(W, sol, pBH, q_s, model)
if isempty(W)||isempty(W(1).lims)
    withinLims = true(numel(sol),1);
else

    nwells     = numel(sol);
    withinLims = true(nwells,1);

    pBH   = double(pBH);
    q_s   = cell2mat( cellfun(@double, q_s, 'UniformOutput', false) );
    q_s = getRatesWOG(q_s, model);

    for wnr = 1:numel(sol)
        lims = W(wnr).lims;
        pBHw  = pBH(wnr);
        q_sw  = q_s(wnr,:);
        qt_sw = sum(q_sw);
        if ~isnumeric(W(wnr).lims)
            if sol(wnr).sign > 0   % injector
                modes   = {'bhp', 'rate'};
                flags = [pBHw > lims.bhp, qt_sw > lims.rate];
            else            % producer
                modes   = {'bhp', 'orat', 'lrat', 'grat', 'wrat', 'vrat'};
                flags = [pBHw           < lims.bhp,  ...
                    q_sw(2)         < lims.orat, ...
                    q_sw(1)+q_sw(2) < lims.lrat, ...
                    q_sw(3)         < lims.grat, ...
                    q_sw(1)         < lims.wrat, ...
                    qt_sw           < lims.vrat];
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

function q_s = getRatesWOG(q_s, model)
switch model
    case 'OW'
        q_s = [q_s, zeros(size(q_s,1), 1)];
    case {'3P', 'BO', 'VO'}
        % nothing
    otherwise
        error(['Unknown model: ', model]);
end
end








