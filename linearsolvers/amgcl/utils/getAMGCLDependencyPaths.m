function [amgclpath, boostpath] = getAMGCLDependencyPaths(varargin)
%Locate AMGCL Build Prerequisites on Current Computer System
%
% SYNOPSIS:
%   [amgclpath, boostpath] = getAMGCLDependencyPaths
%   [amgclpath, boostpath] = getAMGCLDependencyPaths('pn1', pv1, ...)
%
% OPTIONAL PARAMETERS:
%   amgcl_rev - Specific Git revision to download from AMGCL's GitHub
%               repository if necessary.  Character vector.  Default value:
%                  amgcl_rev = 'a551614040f0a7b793b41a4a63386675ca61d8da'.
%
% RETURNS:
%   amgclpath - Location of AMGCL source code.  Character vector.
%
%   boostpath - Location of Boost library headers.  Character vector.
%
% SEE ALSO:
%   `githubDownload`, `amgcl_matlab`, `callAMGCL`.

%{
#COPYRIGHT#
%}

    opt = struct('amgcl_rev', 'a551614040f0a7b793b41a4a63386675ca61d8da');
    opt = merge_options(opt, varargin{:});

    dep_path  = getDependencyFolder();
    amgclpath = get_amgcl_path(dep_path, opt);
    boostpath = get_boost_path(dep_path);

    if ~valid_global_path(amgclpath)
       handle_missing_amgcl_source(amgclpath, dep_path, opt);
    end

    if ~valid_global_path(boostpath)
       handle_missing_boost(boostpath, dep_path)
    end
end

%--------------------------------------------------------------------------

function fldr = getDependencyFolder()
    fldr = fullfile(mrstPath('linearsolvers'), 'amgcl', 'dependencies');
end

%--------------------------------------------------------------------------

function tf = valid_global_path(p)
   tf = ~isempty(p) && ischar(p) && is_directory(p);
end

%--------------------------------------------------------------------------

function dir = get_amgcl_path(dep_path, opt)
   global AMGCLPATH
   dir = AMGCLPATH;

   if isempty(dir)
      dir = fullfile(dep_path, ['amgcl-', opt.amgcl_rev]);
   end
end

%--------------------------------------------------------------------------

function dir = get_boost_path(dep_path)
   global BOOSTPATH
   dir = BOOSTPATH;

   if isempty(dir)
      dir = fullfile(dep_path, 'boost-1_65_1_subset');
   end
end

%--------------------------------------------------------------------------

function handle_missing_amgcl_source(amgclpath, dep_path, opt)
   s = sprintf(['Did not find AMGCL repository path in default ', ...
                'location ''%s'' or in AMGCLPATH global variable. ', ...
                'Would you like to download the files ', ...
                '(approximately 1 MB download)?'], amgclpath);

   if do_download_library(s)
      ensure_prereq_directory_exists(dep_path);

      repo = 'ddemidov/amgcl';
      % Latest tested

      fprintf('Downloading AMGCL...')
      bdir  = fullfile(mrstOutputDirectory(), 'sources', 'github');
      files = githubDownload(repo, 'All', true, 'Base', bdir, ...
                             'Dest', amgclpath, 'Revision', opt.amgcl_rev);

      if isempty(files)
         error('Download:Failure', ...
               'Failed to acquire AMGCL source code from GitHub');
      end

      fprintf(' Ok!\n');
   end
end

%--------------------------------------------------------------------------

function handle_missing_boost(boostpath, dep_path)
   s = sprintf(['Did not find Boost library in default location ', ...
                '''%s'' or in BOOSTPATH global variable. Would you ', ...
                'like to download the the requisite Boost subset ', ...
                '(approximately 1.8 MB download, may take about ', ...
                'one minute to unzip)?'], boostpath);

   if do_download_library(s)
      ensure_prereq_directory_exists(dep_path);

      fprintf('Downloading and extracting Boost. This may take a minute...');

      boost_url = ['https://www.sintef.no/contentassets/', ...
                   '124f261f170947a6bc51dd76aea66129/', ...
                   'boost-1_65_1_subset.zip'];

      unzip(boost_url, dep_path);

      fprintf(' Ok!\n');
   end
end

%--------------------------------------------------------------------------

function ensure_prereq_directory_exists(dname)
   if ~is_directory(dname)
      [ok, msg, id] = mkdir(dname);

      if ~ok
         error(id, ...
              ['Failed to Create AMGCL Prerequisites Directory ', ...
               '''%s'': %s'], dname, msg)
      end
   end
end

%--------------------------------------------------------------------------

function tf = is_directory(elm)
   if exist('isfolder', 'builtin')
      tf = isfolder(elm);
   else
      tf = isdir(elm);
   end
end

%--------------------------------------------------------------------------

function status = do_download_library(msg)
    isDesktop = usejava('desktop');

    if isDesktop
        title = 'Missing dependency';
        choice = questdlg(msg, title, 'Yes', 'No', 'Yes');
    else
        prompt = [msg, ' y/n [y]: '];
        choice = input(prompt, 's');
    end

    status = strcmpi(choice, 'y') || strcmpi(choice, 'yes');
end
