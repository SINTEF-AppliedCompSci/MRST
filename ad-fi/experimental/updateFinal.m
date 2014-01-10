function state = updateFinal(state, state0)
    for i=1:size(state.s,2);
        state.smax(:,i)=max(state.s(:,i),state0.smax(:,i));
        state.smin(:,i)=min(state.s(:,i),state0.smin(:,i));
    end
end