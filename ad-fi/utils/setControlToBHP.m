function schedule = setControlToBHP(schedule)
for cn = 1:numel(schedule.control)
    c = schedule.control(cn);
    for k = 1:size(c.WCONINJE,1)
        schedule.control(cn).WCONINJE{k,4} = 'BHP';
    end
    for k = 1:size(c.WCONPROD,1)
        schedule.control(cn).WCONPROD{k,3} = 'BHP';
    end
end
end
