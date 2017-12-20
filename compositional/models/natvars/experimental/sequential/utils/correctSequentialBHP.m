function states_seq = correctSequentialBHP(modelseq, states_seq, schedule)
    modelseq.transportModel.extraStateOutput = true;
    for i = 1:numel(states_seq)
        i
        W = schedule.control(schedule.step.control(i)).W;
        seq = states_seq{i};
        [~, seq] = transportEquationsNaturalVars(seq, seq, modelseq.transportModel, 1, schedule.control, 'resOnly', true);
        for j = 1:numel(W)
            w = W(j);
            q = seq.wellSol(j).flux;
            c = w.cells;
            WI = w.WI;
            mob = sum(seq.mob(c, :), 2);
            p = seq.pressure(c);

            if ~strcmpi(w.type, 'bhp')
                bhp_new = q./(WI.*mob) - seq.wellSol(j).cdp + p;
                bhp_new = sum(bhp_new.*abs(q))./sum(abs(q));
                states_seq{i}.wellSol(j).bhp = bhp_new;
            end
        end
    end
end