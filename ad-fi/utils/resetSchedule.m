function deck = resetSchedule(deck, smry)
nm = smry.getNms('TIME'); nm = nm{1};
tms = smry.get(nm, 'TIME', :); tms = tms(:);

if isfield(deck.RUNSPEC, 'SI')
   tms = convertFrom(tms, day);
end

tmr = [0; cumsum( deck.SCHEDULE.step.val )];

isRepStep = false(numel(tms), 1);

tmrCount = 1;
for k = 1:numel(tms)
    if tms(k)>tmr(tmrCount)*(1+10*eps)
        isRepStep(k) = true;
        tmrCount = tmrCount +1;
    end
end

repStep = isRepStep(2:end);
deck.SCHEDULE.step.val = diff(tms);
deck.SCHEDULE.step.repStep  = repStep;
deck.SCHEDULE.step.control  = deck.SCHEDULE.step.control(cumsum(repStep));
end
