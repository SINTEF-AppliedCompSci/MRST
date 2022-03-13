function updata = LayeredExactUpscaling(K, krW, krO, Rk, varargin)
% Description
%
% The model consists of three stacked horizontal layers. Bottom and top
% layers are of rock type 1, while the middle layer is of rock type 2.

%{
Copyright 2009-2022 SINTEF Digital, Mathematics & Cybernetics.

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

opt = struct(...
    'Lz',  3, ... % Physical grid dimension in z-direction
    'dz1', 1 ... % Thickness of each of top and bottom layers (are equal)
    );
opt = merge_options(opt, varargin{:});

% Properties

Lz      = opt.Lz;
dz1     = opt.dz1;
perm    = K; % abs perm of rocks, given as [K1, K2]

% Extract properties
dz   = [dz1  (Lz - 2*dz1)];

% Convert data form

sW = krW{1}(:,1);
relperm = nan(numel(sW), 5); % relperm of rocks, [sw krw1 krw2 kro1 kro2]
relperm(:,1) = sW;
relperm(:,2) = krW{1}(:,2);
relperm(:,3) = krW{2}(:,2);
relperm(:,4) = flipud(krO{1}(:,2));
relperm(:,5) = flipud(krO{2}(:,2));

c = Rk{1}(:,1);
rk = nan(numel(c), 3); % rk parameter [cp rk1 rk2]
rk(:,1) = c;
rk(:,2) = Rk{1}(:,2);
rk(:,3) = Rk{2}(:,2);

% Assertions
assert(all(dz) > 0, 'Thickness of layer type 1 too large');
assert(Lz>0 && dz1>0 && all(perm>0));
assert( all(all(relperm>=0)) && all(all(rk>=0)) );



% Upscaling of Absolute Permeability

KU  = upscalePerm(Lz, dz, perm);


% Upscaling of Relative Permeability

krU = upscaleRelperm(Lz, dz, perm, KU, relperm);


% Upscaling of Rk

rkU = upscaleRk(Lz, dz, perm, KU, relperm, krU, rk);


% Set solution structure

% Absperm data
updata.perm = KU;

% Relperm data
updata.krW = cell(1,3);
for d=1:3
    updata.krW{d} = [sW krU.krW(:,d)];
    updata.krO{d} = [sW krU.krO(:,d)]; % TODO: Flip?
end

% Rk data
RkU.val = cell(1,3);
for d=1:3
    RkU.val{d} = rkU.rk(:,:,d);
    RkU.s{d} = rkU.sW;
end
RkU.c = rkU.c;
updata.Rk = RkU;


end


%--------------------------------------------------------------------------
% HELPER FUNCTIONS
%--------------------------------------------------------------------------


function KU = A(Lz, dz, K1, K2)
% Arithmetic average
KU = ( 2*dz(1)*K1 + dz(2)*K2 ) ./ Lz;
end


function KU = H(Lz, dz, K1, K2)
% Harmonic average
KU = Lz ./ ( (2*dz(1))./K1 + dz(2)./K2 );
end


function KU = upscalePerm(Lz, dz, perm)
K1  = perm(1);
K2  = perm(2);
KUxy = A(Lz, dz, K1,K2);
KU   = [KUxy KUxy H(Lz, dz, K1,K2)];
end


function krU = upscaleRelperm(Lz, dz, perm, KU, relperm)

K1  = perm(1);
K2  = perm(2);

krsw = relperm(:,1);
krw1 = relperm(:,2);
krw2 = relperm(:,3);
kro1 = relperm(:,4);
kro2 = relperm(:,5);

sz = [numel(krsw), 3];
krU.sW  = krsw(:);
krU.krW = nan(sz);
krU.krO = nan(sz);

% x- and y-direction (sw* is equal to sw)
krU.krW(:,1) = A(Lz, dz, krw1*K1, krw2*K2) / KU(1);
krU.krW(:,2) = krU.krW(:,1);
krU.krO(:,1) = A(Lz, dz, kro1*K1, kro2*K2) / KU(1);
krU.krO(:,2) = krU.krO(:,1);

% z-direction (need to find sw1 and sw2 first)
sw  = findSaturations(dz, relperm, krsw);
krU.krW(:,3) = H(Lz, dz, interp1(krsw, krw1, sw(:,1))*K1, ...
    interp1(krsw, krw2, sw(:,2))*K2) / KU(3);
krU.krO(:,3) = H(Lz, dz, interp1(krsw, kro1, sw(:,1))*K1, ...
    interp1(krsw, kro2, sw(:,2))*K2) / KU(3);

end


function rkU = upscaleRk(Lz, dz, perm, KU, relperm, krU, rk)

krw1 = relperm(:,2);
krw2 = relperm(:,3);

rkU.sW = relperm(2:end,1); % disregard sw=0
rkU.c  = rk(:,1);
rkU.rk = nan(size(relperm,1)-1, size(rk,1), 3); % nsw x ncp

for i = 1:size(rk,1)
    relperm(:,2) = krw1 ./ rk(i,2);
    relperm(:,3) = krw2 ./ rk(i,3);
    krRkU = upscaleRelperm(Lz, dz, perm, KU, relperm);
    rkUi  = krU.krW ./ krRkU.krW;
    rkU.rk(:, i, :) = rkUi(2:end, :);
end

end


function sw = findSaturations(dz, relperm, swu)
% For each value swu(i), we solve the two equations
%
%   krw1(sw1)*kro2(sw2) = krw2(sw2)*kro1(sw1)
%   2*sw1*dz1 + sw2*dz2 = swu(i)*(2*dz1 + dz2)
%
% to find sw1 and sw2.

krsw = relperm(:,1);
krw1 = relperm(:,2);
krw2 = relperm(:,3);
kro1 = relperm(:,4);
kro2 = relperm(:,5);

sw = nan(numel(swu), 2); % vectors [sw1, sw2]
for i = 1:numel(swu)
    f1 = @(sw1,sw2) interpTable(krsw,krw1,sw1).*...
        interpTable(krsw,kro2,sw2) - ...
        interpTable(krsw,krw2,sw2).*...
        interpTable(krsw,kro1,sw1);
    f2 = @(sw1,sw2) (2*sw1*dz(1) + sw2*dz(2)) - swu(i)*(2*dz(1) + dz(2));
    F  = @(x) [f1(x(1),x(2)); f2(x(1),x(2))];
    x0 = [1;1].*swu(i);
    x  = fsolveadi(F, x0);
    sw(i, :) = x(:)';
end

end



function [x0, converged, iter] = fsolveadi(f, x0)
% Simple fsolve

iterMax = 10;
converged = false;

iter = 0;
while ~converged
    iter = iter + 1;
    x  = initVariablesADI(x0);
    eq = f(x);
    dx = -eq.jac{1}\eq.val;
    x0 = x0 + dx;
    
    if norm(eq.val, Inf) < 1e-12
        converged = true;
    end
    if iter >= iterMax
        disp('**************** MAXIMUM ITERATIONS');
        break
    end
end

end



