function [wellsols,ind]=reMapWellsols(wellsols)
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