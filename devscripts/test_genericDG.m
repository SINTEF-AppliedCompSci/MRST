mrstModule add dg ad-props ad-core ad-blackoil blackoil-sequential ...
    mrst-gui matlab_bgl vemmech vista
mrstVerbose on

%%

gravity reset off
gravity([1,0]);
setup   = getDGTestCase('simple1d', 'n', 3);

setup.modelFV.transportModel.parentModel.fluid.transMult = @(p) 0*p + 1;
setup.modelDG{2}.transportModel.parentModel.fluid.transMult = @(p) 0*p + 1;

%%

% setup.modelFV.transportModel.formulation = 'missingPhase';
[wsFV, stFV, repFV] = simulateScheduleAD(setup.state0, setup.modelFV, setup.schedule);

%%

for dNo = 1:numel(setup.modelDG)
%     setup.modelDG{dNo}.transportModel.formulation = 'missingPhase';
    [wsDG, stDG, repDG] = simulateScheduleAD(setup.state0,setup.modelDG{dNo}, setup.schedule);
end
%%

plotWellSols({wsFV, wsDG})