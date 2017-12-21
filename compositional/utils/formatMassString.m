function s = formatMassString(mass, limit)
% Small utility which returns a human readable string from mass.
    if nargin < 2
        % Second argument can be a limit of time units to output
        limit = inf;
    end
    if sign(mass) > 0
        s = '';
    elseif sign(mass) < 0
        s = '-';
    else
        s = '0 Kilogram';
    end
    mass = abs(mass);

    masssys = {mega*1000*kilogram, 1000*kilogram, kilogram, milli*kilogram};
    massn = {'Megatonne', 'Tonne', 'Kilogram', 'Gram'};
    added = false;
    count = 0;
    for i = 1:numel(masssys)
        a = floor(mass/masssys{i});
        if a > 0
            if added
                space = ', ';
            else
                space = '';
            end
            count = count + added;

            plural = '';
            if a ~= 1, plural = 's'; end
            if count == limit
                % Max depth reached, just print decimal
                s = [s, sprintf([space, '%1.2f ', massn{i}, plural], ...
                                                    mass/masssys{i})]; %#ok
                break
            else
                s = [s, sprintf([space, '%d ', massn{i}, plural], a)]; %#ok
            end

            mass  = mass - a*masssys{i};
            added = true;
        end
    end
end
