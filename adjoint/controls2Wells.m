function [A_N, b_N, A_D, b_D] = controls2Wells(W, schedule, controls)
%
% Create mappings A_N, b_N, A_D, b_D such that
%    q_{tot,N}^n = A_N{n}*u^n + b_N{n}
%    p_{w,D}^n   = A_D{n}*u^n + b_D{n}
% where u^n is the set of controll variables at step n

%{
Copyright 2009-2021 SINTEF Digital, Mathematics & Cybernetics.

This file is part of The MATLAB Reservoir Simulation Toolbox (MRST).

MRST is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

MRST is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with MRST.  If not, see <http://www.gnu.org/licenses/>.
%}


numSteps         = numel(schedule);
numWells         = numel(W);
rateWells     = find( strcmp('rate', {W.type}) );
BHPWells      = find( strcmp('bhp', {W.type}) );
numControls   = numel(controls.well);
% indepControls    = controls.independent;
% numIndepControls = sum(indepControls);
controlWellNums  = [controls.well.wellNum];
nonControlWellNums = setdiff((1:numWells), controlWellNums);

% collect all linEqConst
% ec    = controls.linEqConst;
% numEC = numel(ec);
% ECA = []; ECb = [];
% for k = 1:numEC
%    ECA = [ECA ; ec(k).A];
%    ECb = [ECb ; ec(k).b];
%end

% create mapping u = A u_i + b from independent controls to controls
% indepCont2Cont.A = zeros(numControls, numIndepControls);
% indepCont2Cont.A(indepControls, :) = eye(numIndepControls);
% indepCont2Cont.A(~indepControls,  :) = - ECA(:, ~indepControls)\ECA(:, indepControls);
% indepCont2Cont.b = zeros(numControls, 1);
% indepCont2Cont.b(~indepControls) = ECA(:, ~indepControls)\ECb;

% create mapping w = A u + b from controls to wells
A = accumarray([controlWellNums', (1:numControls)'], ones(numControls,1) , [numWells, numControls]);
b = zeros(numWells, 1);
% Then initial full mapping A u_i + b
% A = cont2Well * indepCont2Cont.A;
% b = cont2Well * indepCont2Cont.b;

% Finally update non-control conditions
for k = 1:numSteps
    vals  =  schedule(k).values;
    b_cur = b;
    b_cur(nonControlWellNums) = vals(nonControlWellNums);

    A_N{k}  = A(rateWells, :);
    b_N{k}  = b_cur(rateWells);

    A_D{k}  = A(BHPWells, :);
    b_D{k}  = b_cur(BHPWells);
end
