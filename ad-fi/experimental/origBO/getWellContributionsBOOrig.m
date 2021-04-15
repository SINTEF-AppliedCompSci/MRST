function [eqs, bWqW, bOqO, bGqG, sol] = getWellContributionsBOOrig(W, pBH, qs, p, rho, b, rs, m, sol)
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

%---------------
maxIt  = 10;
resTol = 0.001;
%---------------


nwell = numel(W);
nPerf = cellfun(@numel, {W.cells})';
perf2well = rldecode((1:nwell)', nPerf);

[eqs{1},eqs{2}, eqs{3}, eqs{4}, bWqW, bOqO, bGqG, alpha, seg_pres] = deal(cell(nwell,1));

for wnr = 1:nwell
    w   = W(wnr);
    s   = sol(wnr);
    pix = perf2well == wnr;

    pBHn = pBH(wnr);
    qWsn = qs{1}(wnr); qOsn = qs{2}(wnr); qGsn = qs{3}(wnr);

    [s, s] = checkLims(s, pBHn, qWsn, qOsn, qGsn);                     %#ok
    res = inf;
    its = 0;

    if false%maxIt > 0
        X = cell(4,1);
        [X{1}, X{2}, X{3}, X{4}] = initVariablesADI(double(pBHn), double(qWsn), double(qOsn), double(qGsn));
        % reservoir props are kept constant in this inner iteration
        %resProps = horzcat(subDb(p, pix), subDb(rho, pix), subDb(b, pix), ...
        %                   subDb(rs, pix), subDb(m, pix));
        resProps = {subDb(p, pix), subDb(rho, pix), subDb(b, pix), ...
                    subDb(rs, pix), subDb(m, pix)};
        while (res > resTol)&&(its < maxIt)
            alpha0 = s.alpha;

            [eqsw, s, s, s, s] = ...
               eqsWellBO(w, X{1}, X(2:4), resProps{:}, s);             %#ok

            dx = SolveEqsADI(eqsw,[]);
            for k = 1:numel(X)
                X{k}.val = X{k}.val + dx{k};
            end
            %X  = cellfun(@(x,y)x.val + y, X, dx);
            res = norm(s.alpha-alpha0, inf);
            [withinLims, s] = checkLims(s, X{:});
            if ~withinLims, res = inf; end
            its = its +1;
        end
        % update
        pBH.val(wnr)   = X{1}.val;
        qs{1}.val(wnr) = X{2}.val;
        qs{2}.val(wnr) = X{3}.val;
        qs{3}.val(wnr) = X{4}.val;
        %update wellSol
        s.pressure = double(X{1});
        s.qWs      = double(X{2});
        s.qOs      = double(X{3});
        s.qGs      = double(X{4});
        %s.sign     = sign(s.qWs+s.qOs+s.qGs);
        sol(wnr) = s;

    end
    sol(wnr) = s;
    %setup equations
    [eqsn, bWqW{wnr}, bOqO{wnr}, bGqG{wnr}, s] = eqsWellBO(w, sub(pBH, wnr), sub(qs, wnr), sub(p,pix), ...
                                        sub(rho, pix), sub(b, pix), sub(rs, pix), sub(m, pix), s);
    for k = 1:4
        eqs{k}{wnr} = eqsn{k};
    end
    alpha{wnr} = s.alpha;
    seg_pres{wnr} = s.seg_pres;
    %
end
% concatenate and set output
for k = 1:4
    eqs{k} = vertcat(eqs{k}{:});
end
bWqW = vertcat(bWqW{:});
bOqO = vertcat(bOqO{:});
bGqG = vertcat(bGqG{:});
end


function sc = sub(c, inx)
if iscell(c)
    sc = cellfun(@(x)x(inx), c, 'UniformOutput', false);
else
    sc = c(inx);
end
end

function sc = subDb(c, inx)
if iscell(c)
    sc = cellfun(@(x)x.val(inx), c, 'UniformOutput', false);
else
    sc = c.val(inx);
end
end




