function [rstrt, rsspec] = readRestartLocal(prefix, kwrds)

try
    rsspec = readEclipseOutputFileUnFmt([prefix, '.RSSPEC']);
catch me
    rsspec = [];
end

if (nargin < 2)
    if ~isempty(rsspec)
        kwrds = unique(rsspec.NAME.values);
        numSteps = numel(rsspec.TIME.values);
    else
        kwrds = {};
        numSteps = 1;
    end
else
    if (nargin == 2)
       numSteps = numel(rsspec.TIME.values);
    else

    end

end

fname = [prefix, '.UNRST'];

[fid, msg] = fopen(fname, 'r', 'ieee-be');
if fid < 0, error([fname, ': ', msg]); end

rstrt = struct();
for k = 1:numel(kwrds)
    rstrt.(genvarname(kwrds{k})) = repmat({{}}, [1, numSteps]);
end

fprintf(['Reading ', num2str(numel(kwrds)), ' fields in ', num2str(numSteps), ' steps\n'])
cnt = 0;
while ~feof(fid)
    [name, field] = readFieldUnFmt(fid);
    if strcmp(name, 'SEQNUM')
        fprintf('*');cnt = cnt+1;
    end
    if ismember(name, kwrds) || isempty(kwrds)
        vnm = genvarname(name);
        if ~isempty(field)
            rstrt.(vnm){cnt} = field.values;
        else
            rstrt.(vnm){cnt} = [];
        end
    end
end
fprintf('\n')
fclose(fid);


