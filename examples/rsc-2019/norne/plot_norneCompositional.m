[ws, states, reports] = getMultiplePackedSimulatorOutputs(problems, 'readFromDisk', false);

%%

ws = cellfun(@(w) {w{1:w.numelData}}, ws, 'unif', false);

%%

loadStates = cell(size(states));
for mNo = 1:2
    st = cell(states{mNo}.numelData,1);
    for sNo = 1:states{mNo}.numelData
        fprintf('Loading state %d of %d ... \n', sNo, states{mNo}.numelData);
        s = states{mNo}{sNo};
        s.transportModel = [];
        st{sNo} = s;
    end
    loadStates{mNo} = st;
end

%%

sd = cellfun(@(s1, s2) compareStates(s1, s2), loadStates{1}(1:19), loadStates{2}, 'unif', false);

%%

plotWellSols(ws)