classdef settingsStruct < handle 
    properties
        allowDL
        dataDirectory
        outputDirectory
        promptDL
        promptMEX
        useMEX
        useOMP
    end
    
    methods
        function val = isfield(settings,fname)
            switch fname
                case 'allowDL'
                    val = true;
                    return
                case 'dataDirectory'
                    val = true;
                    return    
                case 'outputDirectory'
                    val = true;
                    return
                case 'promptDL'
                    val = true;
                    return   
                case 'promptMEX'
                    val = true;
                    return
                case 'useMEX'
                    val = true;
                    return     
                case 'useOMP'
                    val = true;
                    return                     
                otherwise 
                    val = false;
            end
                   
        
        end
    end
    
end
