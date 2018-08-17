function faults = readFaultData()
% Read in fault data from text file
% Each individual fault is loaded as a cell in a cell array

    filename = 'Barents_Sea_faults_CO2_atlas_ED50_UTM32.txt';
    ff = fopen(filename);
    f = 0; % fault index
    while 1
        line = fgetl(ff);
        if line == -1
            fclose(ff);
            break 
        end

        % skip lines that start with # or >
        if line(1) == '#' || line(1) == '>'
            if strcmp(line,'# @D')
                if f > 0
                    faults{f} = [X, Y]; % save previous fault
                end
                f = f+1; % start of new fault
                X = []; Y = [];
            end
            continue
        else
            vals = sscanf(line, '%e');
            assert(numel(vals)==2)
            X = [X; vals(1)];
            Y = [Y; vals(2)];
        end
    end

end