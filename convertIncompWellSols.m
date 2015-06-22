function wellSols = convertIncompWellSols(W, states)
    nstep = numel(states);
    wellSols = cell(nstep, 1);
    
    phases = {'W', 'O'};
    
    for i = 1:nstep
        ws = [];
        for j = 1:numel(W)
            wold = states{i}.wellSol(j);
            
            w.name = W(j).name;
            
            influx = min(wold.flux, 0);
            outflux = max(wold.flux, 0);
            states{i}.s = [states{i}.s 1-states{i}.s];
            for k = 1:2
                fn = ['q', phases{k}, 's'];
                w.(fn) = sum(influx.*states{i}.s(W(j).cells, k)) + sum(outflux*W(j).compi(k));
            end
%             w.qOs = sum(influx.*states(i).s(W(j).cells, 2)) + sum(outflux*W(j).compi(2));
            
            w.bhp = wold.pressure;
            
            ws = [ws; w];
        end
        wellSols{i} = ws;
    end
end
