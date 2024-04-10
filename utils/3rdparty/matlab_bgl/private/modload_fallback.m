function modload
%Fallback Strategy for MATLAB BGL Activation Failure

%{
Copyright 2009-2024 SINTEF Digital, Mathematics & Cybernetics.

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

   s = ['Did not find matlab_bgl library ', ...
                'Would you like to download the files?'];

   if do_download_library(s)
      d = fileparts(fileparts(mfilename('fullpath')));
      run(fullfile(d,'downloadMBGL')) 
   else
       fprintf(2, 'Unable to load matlab_bgl\n');
   end
end



function status = do_download_library(msg)
    status = mrstSettings('get', 'allowDL');
    if status && mrstSettings('get', 'promptDL')
       if mrstPlatform('desktop')
           title = 'Missing dependency';
           choice = questdlg(msg, title, 'Yes', 'No', 'Yes');
       else
           disp(msg);
           choice = input(' y/n [y]: ', 's');
       end
       status = any(strcmpi(choice, {'y', 'yes'}));
   end
end
