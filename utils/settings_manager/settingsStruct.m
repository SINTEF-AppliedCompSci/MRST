classdef settingsStruct < handle & matlab.mixin.SetGetExactNames
%     properties  (Access = private)
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
        
        % Set props
        
        function set.allowDL(obj,value)
            obj.allowDL = value;
        end          

        function set.dataDirectory(obj,value)
            if isdir(value)
                obj.dataDirectory = value;
            else
                error('Set mrstDataDirectory: %s is not a valid directory',value);
            end
        end
        
        function set.outputDirectory(obj,value)
            if isdir(value)
                obj.outputDirectory = value;
            else
                 error('Set mrstOutputDirectory: %s is not a valid directory',value);
            end
        end
        
        
        function set.promptDL(obj,value)
            obj.promptDL = value;
        end
        
        function set.promptMEX(obj,value)
            obj.promptMEX = value;
        end
        
        function set.useMEX(obj,value)
            obj.useMEX = value;
        end
        
        function set.useOMP(obj,value)
            obj.useOMP = value;
        end
        
        function set.useHash(obj,value)
            obj.useHash = value;
        end
        
        % Get props
        
        function val = get.allowDL(obj)
            val = obj.allowDL;
        end
        
        function val = get.dataDirectory(obj)
            val = obj.dataDirectory;
        end
        
        function val = get.outputDirectory(obj)
            val = obj.outputDirectory;
        end        
        
        function val = get.promptDL(obj)
            val = obj.promptDL;
        end
        
        function val = get.promptMEX(obj)
            val = obj.promptMEX;
        end
        
        function val = get.useMEX(obj)
            val = obj.useMEX;
        end
        
        function val = get.useOMP(obj)
            val = obj.useOMP;
        end
        
        function val = get.useHash(obj)
            val = obj.useHash;
        end
        
        
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
         val = settingsStruct.isfield_impl(propertynames(settings), fname);
      end

      function val = isfield_octave(~, fname)
         props = settingsStruct.getSettingsProps();
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

      function [allprops,names,folders] = getSettingsProps()
          % Current list of properties
          names = {'allowDL', 'promptDL', 'promptMEX', 'useMEX', ...
                   'useOMP', 'useHash'};
          folders = {'outputDirectory', 'dataDirectory'};
          allprops = [names,folders];
      end
   end
end

%{
Copyright 2009-2022 SINTEF Digital, Mathematics & Cybernetics.

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
