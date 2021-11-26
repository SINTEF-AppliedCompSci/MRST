function example = makeRandomTraining(example,shutin)

schedule.step  = example.schedule.step;
control        = example.schedule.control;
[ctrlNo,nstep] = rlencode(example.schedule.step.control,1);
inds           = [0; cumsum(nstep)];

ctrlInd = 0;
rng(1);
for m=1:numel(ctrlNo)
    ctrlVals = ceil(0.25*(1:nstep(m))') + ctrlInd;
    schedule.step.control(inds(m)+1:inds(m+1)) = ctrlVals;
    
    % Preserve the first well control
    ctrlInd = ctrlInd+1;
    schedule.control(ctrlInd) = control(ctrlNo(m));

    % Perturb the remaining ones
    for n=ctrlInd+1:max(ctrlVals)
        schedule.control(n) = control(ctrlNo(m));
        W = schedule.control(n).W;
        for i=1:numel(W)
            switch W(i).type
                case 'rate'
                    W(i).val = (.75 + .5*rand)*W(i).val;
                    if ~isempty(W(i).lims) && ~isinf(W(i).lims.rate)
                        W(i).lims.rate = W(i).val;
                    end
                case 'bhp'
                    W(i).val = (.95 + 0.1*rand)*W(i).val;
                    if ~isempty(W(i).lims) && ~isinf(W(i).lims.bhp)
                        W(i).lims.bhp = W(i).val;
                    end
            end
            if shutin && (W(i).sign<0) && (rand(1)>0.9)
                W(i).status = false;
            end
        end
        schedule.control(n).W = W;
    end
    ctrlInd = max(ctrlVals);
end
example.schedule = schedule;
example.name = [example.name '_rand'];