function [wellsols,ind]=reMapWellsols(wellsols)
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

%% assume uniform
ws1=wellsols{1}{1};
ind=nan(numel(ws1),1);
wnames = {ws1.name};
for j=2:numel(wellsols)
    ws2 = wellsols{j}{1};
    wntmp = {ws2.name};
    wnames=intersect(wnames,wntmp);
end

for i=1:numel(wellsols)
    for j=1:numel(wellsols{i})
        ws2=wellsols{i}{j};
        wellsols{i}{j}=ws2(1:numel(wnames));
        wname2={ws2.name};
        for k=1:numel(wnames)
           kk = find([cellfun(@(x) strcmp(wnames{k},x),wname2)]);
           if(~isempty(kk))
                wellsols{i}{j}(k)=ws2(kk);
                ind(k)=kk;
           else
               assert(false);
           end
        end
    end
end
%% clean fields
ff = intersect(fields(wellsols{2}{1}(1)),fields(wellsols{1}{1}(1)));

for i=1:numel(wellsols)
    ws = wellsols{i}    
    for j = 1:numel(ws)
        w = ws{j};
        fn = fields(w);
        rmf = fn(~ismember(fn,ff));
        ws{j} = rmfield(w,rmf);
    end
    wellsols{i} = ws;
end
