function downloadMBGL
   fex = 'http://www.mathworks.com/matlabcentral/fileexchange';
   pth = 'submissions/10922/v/2/download';  % Probably not stable over time
   zip = 'zip';                             % Possibly not stable over time
   url = [fex, '/', pth, '/', zip];

   % Ouput directory.
   dest = fullfile(ROOTDIR, 'utils', '3rdparty', 'matlab_bgl');

   try
      fprintf('Downloading MatlabBGL ZIP archive ... ')
      t = tic;
      unzip(url, dest);
      t = toc(t);
      fprintf('done (%.02f [s])\n\n', t);
   catch
      fallback  = [fex, '/', '10922-matlabbgl'];
      print_url = ['<a href="', fallback, '">', fallback, '</a>'];

      error(['Failed to retrieve MatlabBGL package from known URL\n', ...
             'Please see fall-back location\n  * %s'], print_url);
   end

   if exist('assert', 'builtin'),
      % MatlabBGL contains an 'assert.m' file in its 'test' directory that
      % conflicts with the built-in function of the same name.  Just use
      % the built-in, because the semantics of 'assert.m' are the same as
      % for ASSERT.
      %
      fprintf(['Patching MatlabBGL for compatibility ', ...
               'with recent MATLAB ... ']);
      delete(fullfile(dest, 'matlab_bgl', 'test', 'assert.m'));
      fprintf('done\n');
   end

   % Rewrite the module path for matlab_bgl to point to the correct
   % directory
   fprintf(['Patching ''startup_user'' function for ', ...
            'MatlabBGL installation ... ']);
   [fid, msg] = fopen(fullfile(ROOTDIR, 'startup_user.m'), 'at');
   if fid < 0,
      warning('Open:Fail', 'Failed to open ''startup_user.m'': %s', msg);
   end
   fprintf(fid, 'mrstPath reregister matlab_bgl ...\n    %s\n', ...
           fullfile(dest, 'matlab_bgl'));
   fclose(fid);
   fprintf('done\n');

   % Create a modload.m that's compatible with the new, third-party
   % component.
   copyfile(fullfile(dest, 'private', 'modload.m.in'), ...
            fullfile(dest, 'private', 'modload.m'));

   clear modload modload_fallback
end
