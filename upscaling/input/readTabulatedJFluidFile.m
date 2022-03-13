function T = readTabulatedJFluidFile(masterfluid)
% Read tabulated J fluid. See examples for syntax.

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


% Read masterfluid file

[fid, msg] = fopen(masterfluid);
if fid < 0,
   error(msgid('Open:Failure'), ...
         'Failed to open rock-list file ''%s'': %s.', masterfluid, msg);
end

% read number of regions
numfiles = str2double(fgetl(fid));


done = false;
   data = [];
   i = 0;
   while ~done, % test on numfiles also?
      i = i + 1;
      lin  = fgetl(fid);
      lin  = regexprep(lin, '--.*$', '');  % Remove comments/disabled data.
      lin  = regexprep(lin, '#.*$', '');  % Remove comments/disabled data.
       lin  = regexprep(lin, '%.*$', '');  % Remove comments/disabled data.
      done = feof(fid) || i>=numfiles ;              % Record complete?
      data = [data, ' ', lin];  %#ok       % Willfully ignore MLINT.
   end

fluidlist = regexp(strtrim(data), '\s+', 'split');

fclose(fid);

% Read fluid files

T = cell([numfiles, 1]);
fluiddir = fileparts(masterfluid);

if numel(fluidlist) == numfiles   % scalar relperm format

   for t = 1 : numfiles,
      fn = fullfile(fluiddir, fluidlist{t});
      [fid, msg] = fopen(fn);
      if fid < 0,
         warning(msgid('Open:Failure'), ...
                ['Failed to open tabulated saturation function ', ...
                 '''%s'': %s. Continuing...'], fn, msg);
         continue
      end
      % (fid, ntab, ncol)
      t_w = readRelPermTableStatoil(fid, 4);
      % HACK: set first value in table to be connate water saturation, krw = 0
      t_w(1,2) = 0;
      t_w(end,3) = 0;

      T{t} = t_w;
      fclose(fid);
   end

elseif numel(fluidlist) == 2*numfiles   % tensor relperm format - use krx
   T = cell([numfiles, 1]);
   for i = 1:2:numel(fluidlist)
      % Read water relperm table
      fn = fullfile(fluiddir, fluidlist{i});
      [fid, msg] = fopen(fn);
      if fid < 0,
         warning(msgid('Open:Failure'), ...
                ['Failed to open tabulated saturation function ', ...
                 '''%s'': %s. Continuing...'], fn, msg);
         continue
      end
      % (fid, ntab, ncol)
      t_w = readRelPermTableStatoil(fid, 5);
      fclose(fid);
      % HACK: set first value in table to be connate water saturation, krw = 0
      t_w(1,3) = 0;

      % Read oil relperm table
      fid = fopen(fullfile(fluiddir, fluidlist{i+1})) ;
      t_o = readRelPermTableStatoil(fid, 5);
      % HACK: set first value in table to kro = 0, oil immobile at max
      % watersat
      t_o(end,3) = 0;
      fclose(fid);
      % s_w      s_w      krw     kro          J_pc
      T{(i+1)/2} = [t_w(:,2) t_w(:,3) t_o(:,3) t_w(:,1)];
   end

else
   disp('FAILED')
end

end


%--------------------------------------------------------------------------
% Helpers follow.
%--------------------------------------------------------------------------
function T = readRelPermTableStatoil(fid, ncol)
      % A rel-perm table is textually represented as a sequence of rows,
      % each row containing 'ncol' numbers (or defaulted elements of the
      % form '1*').  The sequence is terminated by a single slash ('/')
      % character.
      %
      % Read data into a character string, split on whitespace, convert to
      % DOUBLE, and (if required) insert default values computed by linear
      % interpolation.
      %
      data = readRecordStringTabulated(fid);    % Treat table data as single record.
      data = regexp(strtrim(data), '\s+', 'split');

      T = convertTable(reshape(data, ncol, []) .');
end

function data  = readRecordStringTabulated(fid)
done = false;
   data = [];

   while ~done,
      lin  = fgetl(fid);
      if lin ~= -1
      lin  = regexprep(lin, '--.*$', '');  % Remove comments/disabled data.
      lin  = regexprep(lin, '#.*$', '');  % Remove comments/disabled data.
       lin  = regexprep(lin, '%.*$', '');
      done = any(lin == '/') || feof(fid) ;              % Record complete?
      data = [data, ' ', lin];  %#ok       % Willfully ignore MLINT.
      else
         done = true;
      end
   end
 %  data = data(1 : find(data == '/') - 1); % Don't include terminator char.
  % data(data == '''') = '';                % Exclude any quote characters.
end

function T = convertTable(T)
   T = cellfun(@(s) sscanf(s, '%f'), strrep(T, '1*', 'NaN'));
end
