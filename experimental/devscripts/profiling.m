if 0
    [wsDG, stDG, repDG] = sim('modelDG', 2);
else
    ix = 1:10;
    schdl = setup.schedule;
    schdl.step.val = schdl.step.val(ix);
    schdl.step.control = schdl.step.control(ix);
    [wsDG, stDG, repDG] = simulateScheduleAD(setup.state0, setup.modelDG{2}, schdl);
end