function val = computeStatistic(vals, stat, prop)
type = 1; % standard mean/var
if nargin == 3
    if strcmp('TOF', prop(1:3))
        type = 2; % harmonic mean/var for TOFs
    end
end

n = size(vals,2);
switch stat
    case 'mean'
        if type == 1
            val = mean(vals, 2);
        else % harmonic
            val = 1./(mean(1./vals, 2));
        end
    case 'std'
        if type == 1
            val = std(vals, 0, 2);
        else % harmonic
            val = 1./(std(1./vals, 0, 2)); % or something ...
        end
    case 'max diff'
        val = max(vals, [], 2) - min(vals, [], 2);
end
end
        
        