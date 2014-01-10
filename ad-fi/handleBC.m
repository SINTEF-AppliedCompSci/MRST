function eq = handleBC(W, pBHP, qWs, qOs, qGs, scalFacs)
if nargin < 6
    scalFacs.rate = 1; scalFacs.pressure = 1;
end

eq = pBHP;
isInj = (vertcat(W(:).sign)>0);
compi = vertcat(W(:).compi);

% bhp (injector or producer)
inx = find(cellfun(@(x)strcmp('bhp',x), {W.type}'));
if ~isempty(inx)
    val = vertcat(W(inx).val);
    eq(inx) = (pBHP(inx) - val)/scalFacs.pressure;
end

%rate (injector)
inx = find(cellfun(@(x)strcmp('rate',x), {W.type}'));
if ~isempty(inx)
    assert(all(isInj(inx)));
    inxW = inx(logical(compi(inx,1)));
    inxO = inx(logical(compi(inx,2)));
    inxG = inx(logical(compi(inx,3)));
    if ~isempty(inxW)
        val = vertcat(W(inxW).val);
        eq(inxW) = (qWs(inxW)-val)/scalFacs.rate;
    end
    if ~isempty(inxG)
        val = vertcat(W(inxG).val);
        eq(inxG) = (qGs(inxG)-val)/scalFacs.rate;
    end
    if ~isempty(inxO)
        if (all(val == 0)) && strcmp(W.type, 'rate')
            % this is a dummy well!
        else
            warning('handleBC:oil_inj','Current setup will lead to injection of oil.');
        end
        val = vertcat(W(inxO).val);
        eq(inxO) = (qOs(inxO)-val)/scalFacs.rate;
    end
end

%orat (producer)
inx = find(cellfun(@(x)strcmp('orat',x), {W.type}'));
if ~isempty(inx)
    assert(all(~isInj(inx)));
    val = vertcat(W(inx).val);
    eq(inx) = (qOs(inx)-val)/scalFacs.rate;
end

%wrat (producer)
inx = find(cellfun(@(x)strcmp('wrat',x), {W.type}'));
if ~isempty(inx)
    assert(all(~isInj(inx)));
    val = vertcat(W(inx).val);
    eq(inx) = (qWs(inx)-val)/scalFacs.rate;
end

%grat (producer)
inx = find(cellfun(@(x)strcmp('grat',x), {W.type}'));
if ~isempty(inx)
    assert(all(~isInj(inx)));
    val = vertcat(W(inx).val);
    eq(inx) = (qGs(inx)-val)/scalFacs.rate;
end

%lrat (producer)
inx = find(cellfun(@(x)strcmp('lrat',x), {W.type}'));
if ~isempty(inx)
    assert(all(~isInj(inx)));
    val = vertcat(W(inx).val);
    eq(inx) = (qOs(inx)+qWs(inx)-val)/scalFacs.rate;
end

end
