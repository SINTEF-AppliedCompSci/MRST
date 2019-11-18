function aq = createAquiSol(rstrt, repStep, ijk2act)
[aq, unitNo]   = getRestartAquiInfo(rstrt, repStep);

if ~isempty(aq)
    unit = {'metric', 'field', 'lab'};
    unit = unit{unitNo};
    u    = getUnitSystem(unit);
    
    for k = 1:numel(aq)
        aq(k).cells = ijk2act(aq(k).cijk);
        aq(k).p     = convertFrom(aq(k).p, u.press);
        aq(k).qW    = convertFrom(aq(k).qW,  u.resvolume/u.time);
        aq(k).flux  = convertFrom(aq(k).flux, u.resvolume/u.time);
    end
    
    aq = rmfield(aq, 'cijk');
end
end
