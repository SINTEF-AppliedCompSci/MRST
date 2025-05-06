function IFT = eCPAInterfacialTension(eos, T, x, IFT_pure,L)
% Estimate interfacial tension of mixtures

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

R = 8.314462618;
names = eos.ECPACompositionalMixture.names;
nmole = eos.getNumberOfMolecules;
[M, I] = getIFTparameters(names);

c = 0;IFT = 0;
L = value(L);
singlePhase = L==1 | L==0;
if iscell(x)
    IFT_pure = expandMatrixToCell(IFT_pure);
    for i = 1 : nmole
        for j = 1 : nmole
            if i ~= j
                c = c + x{i}.*x{j}.*M(i,j).*(298.15./T).^(I(i,j)-1);
            end
        end
        IFT = IFT + x{i} .* IFT_pure{i};
    end
    IFT = IFT - 1e-3 .* R .* T .* c;
    IFT(singlePhase) = IFT(singlePhase).*0;
else
    for i = 1 : nmole
        for j = 1 : nmole
            if i ~= j
                c = c + x(:,i).*x(:,j).*M(i,j).*(298.15./T).^(I(i,j)-1);
            end
        end
    end
    IFT = sum(x(:,1:nmole).*IFT_pure, 2) - 1e-3 .* R .* T .* c;
    IFT(singlePhase) = 0;
end

end

function [M, I] = getIFTparameters(names)
ncomp = numel(names);
[M, I] = deal(zeros(ncomp,ncomp));
for i = 1 : ncomp - 1
    switch(lower(names{i}))
        case {'water', 'h2o'}
            for j = i + 1 : ncomp
                switch(lower(names{j}))
                    case {'carbondioxide', 'co2'}
                        M(i,j) = 0.356; M(j,i) = M(i,j);
                        I(i,j) = 3.049; I(j,i) = I(i,j);
                    case {'hydrogensulfide', 'h2s'}
                        M(i,j) = 0.323; M(j,i) = M(i,j);
                        I(i,j) = 4.666;  I(j,i) = I(i,j);
                    case {'methane', 'ch4'}
                        M(i,j) = 1.534; M(j,i) = M(i,j);
                        I(i,j) = 4.046;  I(j,i) = I(i,j);
                end
            end
        case {'carbondioxide', 'co2'}
            for j = i + 1 : ncomp
                switch(lower(names{j}))
                    case {'methane', 'ch4'}
                        M(i,j) = 0.006; M(j,i) = M(i,j);
                        I(i,j) = 3.799;  I(j,i) = I(i,j);
                end
            end
    end
end

end
