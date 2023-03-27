function state = updateFinalVE(state, state0,f,G)
    for i=1:size(state.s,2);
        state.smax(:,i)=max(state.s(:,i),state0.smax(:,i));
        state.smin(:,i)=min(state.s(:,i),state0.smin(:,i));  
    end
    if isfield(f,'dis_rate')
        min_rs= minRs(state.pressure,state.s(:,2),state.sGmax,f,G);
        min_rs=min_rs./state.s(:,1);
        assert(all(state.rs - min_rs >= -2.0e-6)); % @@
    end
          
end