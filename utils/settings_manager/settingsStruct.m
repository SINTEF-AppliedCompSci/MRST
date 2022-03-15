classdef settingsStruct < handle
    properties
        allowDL
        dataDirectory
        outputDirectory
        promptDL
        promptMEX
        useMEX
        useOMP
        useHash
    end

    methods
        function val = isfield(settings, fname)
           if mrstPlatform('octave')
              val = settings.isfield_octave(fname);
           else
              val = settings.isfield_matlab(fname);
           end
        end
    end

   methods (Access = private)

      function val = isfield_matlab(settings, fname)
         val = settingsStruct.isfield_impl(properties(settings), fname);
      end

      function val = isfield_octave(~, fname)

         props = getSettingsProps(settings);
         val = settingsStruct.isfield_impl(props, fname);
      end
          
   end

   methods (Access = private, Static)
      function val = isfield_impl(fields, fname)
         fname = cellstr(fname);
         [i, j] = blockDiagIndex(numel(fname), numel(fields));

         m = strcmp(reshape(fname(i), numel(fname), []), ...
                    reshape(fields(j), numel(fname), []));

         val = reshape(any(m, 2), size(fname));
      end
      
             
      function [allprops,names,folders] = getSettingsProps(settings)
          % Current list of properties
           names = {'allowDL', 'promptDL', 'promptMEX', 'useMEX', ...
               'useOMP', 'useHash'};
          folders = {'outputDirectory', 'dataDirectory'};
          allprops = [names,folders];
      end      
      
      
   end
end
