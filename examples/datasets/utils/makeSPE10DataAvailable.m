function ok = makeSPE10DataAvailable()
%Ensure availability of Models 1 and 2 from tenth SPE Comparative Solution Project
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

   ok =        output_exists(model1_data()) || download_model1();
   ok = ok && (output_exists(model2_matfile_name()) || download_model2());
end

%--------------------------------------------------------------------------

function ok = download_model1()
   assert (~ output_exists(model1_data()), 'Internal Error');

   ok = have_model1_perm_kr() || do_download_model1();

   if ok
      % We should *always* end up here.

      ok = write_model1_matfile();
   end
end

%--------------------------------------------------------------------------

function ok = download_model2()
   assert (~ output_exists(model2_matfile_name()), 'Internal Error');

   ok = have_perm_poro_input() || do_download_model2();

   if ok
      % We should *always* end up here.

      ok = write_model2_rock_matfile();
   end
end

%--------------------------------------------------------------------------

function files = model1_data()
   files = { model1_matfile_name(), 'krog.txt' };
end

%--------------------------------------------------------------------------

function fname = model1_matfile_name()
   fname = 'model1_data.mat';
end

%--------------------------------------------------------------------------

function fname = model2_matfile_name()
   fname = 'spe10_rock.mat';
end

%--------------------------------------------------------------------------

function ok = do_download_model1()
   % Data files not available in unpacked form.  Extract local
   % archive *or* the archive file downloaded from the SPE web site and
   % retrieve relative permeability data.

   do_download('perm_case1.zip', 'first');

   wget = mrstWebSave();

   sgof = wget(output_file('krog.txt'), ...
               [csp_url(), 'rel_perm_tab.txt']);

   ok = ~isempty(sgof) && have_model1_perm_kr();
end

%--------------------------------------------------------------------------

function ok = do_download_model2()
   % Data files not available in unpacked form.  Extract local
   % archive *or* the archive file downloaded from the SPE web site.

   do_download('por_perm_case2a.zip', 'second');

   ok = have_perm_poro_input();
end

%--------------------------------------------------------------------------

function do_download(zipfile, comment)
   if output_exists(zipfile)
      url = output_file(zipfile);
   else
      url = [csp_url(), zipfile];
   end

   dispif(mrstVerbose, ...
         ['Please wait while the %s SPE10 ', ...
          'dataset is downloaded...'], comment);

   unzip(url, output_directory());
   dispif(mrstVerbose, 'Done\n');
end

%--------------------------------------------------------------------------

function ok = write_model1_matfile()
   assert (have_model1_perm_kr(), 'Internal Error');

   % Heterogeneous permeability field, homogeneous porosity field.
   perm = load_model_data('perm_case1.dat', perm_filter(), 3);
   poro = repmat(0.2, [ size(perm, 1), 1 ]);

   rock = struct('perm', perm, 'poro', poro);

   ncell = 100 * 1 * 20;

   ok = all([size(rock.perm, 1), size(rock.poro, 1)] == ncell);

   kr_deck = model1_deck();                                     %#ok<NASGU>

   save(output_file(model1_matfile_name()), 'rock', 'kr_deck');
end

%--------------------------------------------------------------------------

function ok = write_model2_rock_matfile()
   assert (have_perm_poro_input(), 'Internal Error');

   % Permeability data.
   rock.perm = load_model2_data('perm', perm_filter(), 3);

   % Porosity data.
   rock.poro = load_model2_data('phi', identity_map(), 1);

   % Verify size of input data.
   %
   ncell = 60 * 220 * 85;

   ok = all([size(rock.perm, 1), numel(rock.poro)] == ncell);

   % Save data in form more amenable to subsequent M processing.
   %
   save(output_file(model2_matfile_name()), 'rock')
end

%--------------------------------------------------------------------------

function deck = model1_deck()
%Stripped-down simulation deck for ADI-fluid (relperm) construction
%
% This input-deck contains just enough information to construct relative
% permeability curves from the benchmark data.  The 'krog.txt' is formatted
% like ECLIPSE's SGOF keyword (first three lines provide human context).

   [fid, msg] = fopen(output_file('krog.txt'), 'rt');

   if fid < 0
      error('Fopen:Fail', 'Failed to Open O/G Rel-Perm File: %s', msg);
   end

   % TEXTSCAN available since R14 (MATLAB 7.0).
   krdata = textscan(fid, '%f %f %f %f', ...
                     'HeaderLines'  , 3, ...
                     'CollectOutput', true);

   fclose(fid);

   deck = struct('GRID'   , struct()                  , ...
                 'PROPS'  , struct('SGOF', { krdata }), ...
                 'REGIONS', struct());
end

%--------------------------------------------------------------------------

function x = load_model2_data(quantity, filter, ncol)
   x = load_model_data(['spe_', quantity, '.dat'], filter, ncol);
end

%--------------------------------------------------------------------------

function x = load_model_data(file, filter, ncol)
   [fid, msg] = fopen(output_file(file), 'rt');

   if fid < 0, error(msg); end

   x = reshape(filter(fscanf(fid, '%f')), [], ncol);

   fclose(fid);
end

%--------------------------------------------------------------------------

function f = perm_filter()
   f = @(perm) convertFrom(perm, milli*darcy);
end

%--------------------------------------------------------------------------

function i = identity_map()
   i = @(x) x;
end

%--------------------------------------------------------------------------

function tf = have_model1_perm_kr()
   tf = output_exists({ 'perm_case1.dat', 'krog.txt' });
end

%--------------------------------------------------------------------------

function tf = have_perm_poro_input()
   tf = output_exists(strcat('spe_', {'perm', 'phi'}, '.dat'));
end

%--------------------------------------------------------------------------

function tf = output_exists(fname)
   if ischar(fname), fname = { fname }; end

   tf = all(cellfun(@(fn) exist(fn, 'file') == 2, output_file(fname)));
end

%--------------------------------------------------------------------------

function fname = output_file(fname)
   fname = fullfile(output_directory(), fname);
end

%--------------------------------------------------------------------------

function odir = output_directory()
   odir = getDatasetPath('spe10', 'skipAvailableCheck', true);
end

%--------------------------------------------------------------------------

function url = csp_url()
   url = 'http://www.spe.org/web/csp/datasets/';
end
