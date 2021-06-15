function eqs = addMinearals(eqs, I, M, T, state, state0, G,  W, fluid, fbG, s, dt, dpG, mobG, bG, sG, sG0, bGqG, pvMult0, pvMult)
%Undocumented Utility Function

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

eqn_num=numel(eqs);
WIinj=vertcat(W.I);
eqn_numI=eqn_num;
wc=[W.cells]';
assert(size(fluid.ILn,2)== size(state.I,2));
assert(size(fluid.IRn,2)== size(state.I,2));
for i=1:size(state.I,2)
    eqn_num=eqn_numI+i;
    bIqI = I{i}(wc).*bGqG;
    ind=bGqG<0;
    if(any(ind))
        bIqI(ind)  = WIinj(ind,i).*bGqG;
    end
    %upc = (double(dpG)>=0);
    upc = (double(dpG)>=0);
    bIvI = s.faceUpstr(upc, I{i}.*bG.*mobG).*s.T.*dpG;
    eqs{eqn_num } = (s.pv/dt).*( pvMult.*bG.*sG.*I{i} - pvMult0.*fbG(state0.pressure).*sG0.*state0.I(:,i)) + s.div(bIvI);
    eqs{eqn_num}(wc) = eqs{eqn_num}(wc) + bIqI;
end
eqn_numM=eqn_num;
assert(size(state.M,2)==size(fluid.Mn,2));
for i=1:size(state.M,2)
    eqn_num=eqn_numM+i;
    eqs{eqn_num } =(s.pv/dt).*(M{i} - state0.M(:,i));
end

fluidL_rate=cell(size(fluid.IRn,1),1);
fluidR_rate=cell(size(fluid.IRn,1),1);
for j=1:size(fluid.IRn,1)
   fluidL_rate{j}=fluid.LR{j}(T).*sG;
   fluidR_rate{j}=fluid.RR{j}(T).*sG;
end


% ADD REACTIONS
for j=1:size(fluid.IRn,1)
    % reaction j
    lograteR=zeros(G.cells.num,1);
    lograteL=zeros(G.cells.num,1);
    % add eps to avoid log of zero
    for i=1:size(fluid.IRn,2)
        lograteR=lograteR+fluid.IRn(j,i).*log(I{i}+eps);
    end
    for i=1:size(fluid.ILn,2)
        lograteL=lograteL+fluid.ILn(j,i).*log(I{i}+eps);
    end
    for i=1:size(fluid.Mn,2)
        lograteR=lograteR+fluid.Mn(j,i).*log(M{i}+eps);
    end
    for i=1:size(fluid.ILn,2)
       %if( fluid.ILn(j,i)>0)% handle righ  side of equations
       eqs{eqn_numI+i}=eqs{eqn_numI+i}-fluidL_rate{j}.*fluid.IRn(j,i).*exp(lograteL); % to the right
       eqs{eqn_numI+i}=eqs{eqn_numI+i}+fluidR_rate{j}.*fluid.IRn(j,i).*exp(lograteR); % from the right
       %end
    end
    for i=1:size(fluid.IRn,2) %handle  left side of equation
       eqs{eqn_numI+i}=eqs{eqn_numI+i}-fluidR_rate{j}.*fluid.ILn(j,i).*exp(lograteR); % to the left
       eqs{eqn_numI+i}=eqs{eqn_numI+i}+fluidL_rate{j}.*fluid.ILn(j,i).*exp(lograteL); % from the left
    end
    for i=1:size(fluid.Mn,2)
       eqs{eqn_numM+i}=eqs{eqn_numM+i}-fluidL_rate{j}.*fluid.Mn(j,i).*exp(lograteL);
       eqs{eqn_numM+i}=eqs{eqn_numM+i}+fluidR_rate{j}.*fluid.Mn(j,i).*exp(lograteR);
    end
end

end
