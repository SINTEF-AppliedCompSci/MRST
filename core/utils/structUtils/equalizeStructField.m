function [mstruct, bothUnAssigned] = equalizeStructField(mstruct, fieldnamelist1, fieldnamelist2, varargin)

    value1 = getStructField(mstruct, fieldnamelist1);
    value2 = getStructField(mstruct, fieldnamelist2);

    bothUnAssigned = false;
    
    if isUnAssigned(value1)

        if isUnAssigned(value2)
            
            bothUnAssigned = true;
            return
            
        else
            
            mstruct = setStructField(mstruct, fieldnamelist1, value2);

        end

    else

        if isUnAssigned(value2)
            
            mstruct = setStructField(mstruct, fieldnamelist2, value1);

        else

            if isequal(value1, value2)

                return

            else

                opt = struct('force', 'false', ...
                             'warn', true);
                opt = merge_options(opt, varargin{:});

                if opt.force
                    errorMessage = sprintf('Different values are given for the fields mstruct.%s and mstruct.%s. We do not know which one to choose...', ...
                                           getPrintableName(fieldnamelist1)                                                                                  , ...
                                           getPrintableName(fieldnamelist2));
                    mstruct = setStructField(mstruct, fieldnamelist2, value1, 'handleMisMatch', 'error', 'errorMessage', errorMessage);
                    if opt.warn
                        fprintf('Fist value given in equalizeStructField is taken\n');
                    end
                else
                    error('mismatch in equalizeStructField');
                end
                
            end
        end

    end
      
end


function namestr = getPrintableName(fieldnamelist)

    if ischar(fieldnamelist)

        namestr = getPrintableName({fieldnamelist})
        return
        
    end

    namestr = strjoin(fieldnamelist, '.')

end
