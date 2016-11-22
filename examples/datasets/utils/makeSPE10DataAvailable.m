function ok = makeSPE10DataAvailable()
%Ensure availability of Model 2 from tenth SPE Comparative Solutions Project
%
% SYNOPSIS:
%   ok = makeSPE10DataAvailable
%
% DESCRIPTION:
%   If the dataset is not already present on disk in the directory named by
%
%       getDatasetPath('spe10')
%
%   then this function will attempt to download the petrophysical
%   properties (permeability and porosity) from the SPE website into this
%   directory.
%
% PARAMETERS:
%   None.  This function operates on constant data.
%
% RETURNS:
%   ok - Whether or not the petrophysical properties of Model 2 from the
%        tenth SPE Comparative Solutions Project is availble in the SPE 10
%        dataset directory.
%
% NOTE:
%   The on-disk representation is a 'rock' structure stored in the file
%   'spe10_rock.mat' in the directory identified by
%
%       getDatasetPath('spe10')
%
%   The permeability data is stored in SI units (metres squared).
%
% SEE ALSO:
%   getDatasetPath.

%{
#COPYRIGHT#
%}

   ok = output_file_exists(matfile_name()) || download();
end

%--------------------------------------------------------------------------

function ok = download()
   assert (~ output_file_exists(matfile_name()), 'Internal Error');

   ok = have_perm_poro_input() || do_download();

   if ok,
      % We should *always* end up here.

      ok = write_rock_matfile();
   end
end

%--------------------------------------------------------------------------

function fname = matfile_name()
   fname = 'spe10_rock.mat';
end

%--------------------------------------------------------------------------

function ok = do_download()
   % Data files not available in unpacked form.  Extract local
   % archive *or* the archive file downloaded from the SPE web site.

   zipfile = 'por_perm_case2a.zip';
   if output_file_exists(zipfile),
      url = output_file(zipfile);
   else
      url = ['http://www.spe.org/web/csp/datasets/', zipfile];
   end

   dispif(mrstVerbose, ...
         ['Please wait while the second SPE10 ', ...
          'dataset is downloaded...']);

   unzip(url, output_directory());
   dispif(mrstVerbose, 'Done\n');

   ok = have_perm_poro_input();
end

%--------------------------------------------------------------------------

function ok = write_rock_matfile()
   assert (have_perm_poro_input(), 'Internal Error');

   % Permeability data.
   rock.perm = load_data('perm', @(k) convertFrom(k, milli*darcy), 3);

   % Porosity data.
   rock.poro = load_data('phi', @(phi) phi, 1);

   % Verify size of input data.
   %
   ncell = 60 * 220 * 85;

   ok = all([size(rock.perm, 1), numel(rock.poro)] == ncell);

   % Save data in form more amenable to subsequent M processing.
   %
   save(output_file(matfile_name()), 'rock')
end

%--------------------------------------------------------------------------

function x = load_data(quantity, filter, ncol)
   [fid, msg] = fopen(output_file(['spe_', quantity, '.dat']), 'rt');

   if fid < 0, error(msg); end

   x = reshape(filter(fscanf(fid, '%f')), [], ncol);

   fclose(fid);
end

%--------------------------------------------------------------------------

function tf = have_perm_poro_input()
   tf = output_file_exists('spe_perm.dat') && ...
        output_file_exists('spe_phi.dat');
end

%--------------------------------------------------------------------------

function tf = output_file_exists(fname)
   tf = exist(output_file(fname), 'file') == 2;
end

%--------------------------------------------------------------------------

function fname = output_file(fname)
   fname = fullfile(output_directory(), fname);
end

%--------------------------------------------------------------------------

function odir = output_directory()
   odir = getDatasetPath('spe10', 'skipAvailableCheck', true);
end
