function schedule = scheduleFromSummary(schedule, summary)
    % There are n-1 control steps for n summary steps
    n_t = numel(schedule.step.control);
    schedule.control = repmat(schedule.control(1), n_t, 1);
    schedule.step.control = (1:n_t) .';

    for t = 1:n_t
        wconprod = schedule.control(t).WCONPROD;
        for w = 1:numel(wconprod(:, 1))
            name = wconprod(w,1);
            assert(strcmpi(wconprod(w, 3), 'bhp'));
            schedule.control(t).WCONPROD{w,9} = convertFrom(summary.get(name, 'WBHP', t+1), psia);
        end

        if 1
        wconinje = schedule.control(t).WCONINJE;
        for w = 1:numel(wconinje(:, 1))
            name = wconinje(w,1);
            assert(strcmpi(wconinje(w, 4), 'bhp'));

            schedule.control(t).WCONINJE{w,7} = convertFrom(summary.get(name, 'WBHP', t+1), psia);

            schedule.control(t).WCONINJE{w,5} = inf;
        end
        end
    end

end
