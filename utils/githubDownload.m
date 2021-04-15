function files = githubDownload(repository, varargin)
%Download objects from GitHub (.ZIP or collection of files)
%
% SYNOPSIS:
%   files = githubDownload(repository, 'pn1', pv1, ...)
%
% PARAMETERS:
%   repository - Name of repository from which to construct object URLs.
%
% OPTIONAL PARAMETERS:
%  'Revision' - Git revision of objects.  String.  Subject to expansion
%               through the equivalent of 'rev-parse'.
%               Default value: Revision = 'master'.
%
%  'File' -     Name of file or files to download from GitHub. String or
%               cell array of strings respectively.
%               Default value: File = {} (No specific filename).
%
%  'All' -      Whether or not to download the entire contents of the
%               GitHub repository (at specified revision).
%               Logical.  Default value All=false.
%
%  'Base' -     Base directory in which a subdirectory named after the
%               repository will be created to host the requisite contents
%               of the GitHub repository.
%               String.  Default value: Base = mrstDataDirectory().
%
%  'Dest' -     Directory into which the downloaded file set will be moved.
%               Only taken into account if non-empty. 
%               String.  Default value: Dest = '' (leave downloaded files
%               in original download location).
%
%  'Pause' -    Amount of time to wait between successive requests to
%               GitHub web services.  Only relevant when explicit file list
%               has more than one entry.
%               Non-negative scalar.  Default value: Pause = 5 sec.
%
% NOTES:
%   The user must specify either an explicit list of files or option 'All'.
%
%   If option 'All' is set, then this takes precedence over any explicit
%   file list.  In other words, an explicit file list is ignored if option
%   'All' is set.
%
%   Using a small value for 'Pause' increases the likelihood that the next
%   request will fail.  Pause should usually be at least two seconds.
%
% RETURNS:
%   files - Cell array of strings containing file names, specific to local
%           computer system, of the objects downloaded from GitHub.
%
% EXAMPLE:
%   repo     = 'OPM/opm-data';
%   rev      = '2198d5b';
%   filebase = 'flow_diagnostic_test/eclipse-simulation';
%   files    = strcat('SIMPLE_SUMMARY.', { 'INIT', 'UNRST', 'UNSMRY' });
%
%   files = githubDownload(repo, 'Revision', rev, ...
%                          'File', strcat(filebase, '/', files));
%
% SEE ALSO:
%   `websave`, `unzip`, `mrstDataDirectory`.

%{
Copyright 2009-2021 SINTEF Digital, Mathematics & Cybernetics.

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

   opt = struct('revision', 'master',          ...
                'base',     mrstDataDirectory, ...
                'dest',     '',                ...
                'file',     {{}},              ...
                'pause',    default_wait,      ...
                'all',      false);

   opt = merge_options(opt, varargin{:});

   assert (opt.all || ~isempty(opt.file), ...
           'Must either specify file list or all files');

   if opt.all
      download = @download_zip;
      move     = @(files) move_all_to_dest(files, repository, opt);
   else
      download = @download_files;
      move     = @(files) move_fileset_to_dest(files, opt);
   end

   files = download(repository, object_url(repository, opt), opt);

   if ~ isempty(opt.dest)
      files = move(files);
   end
end

%--------------------------------------------------------------------------

function wait = default_wait
   wait = 5*second;
end

%--------------------------------------------------------------------------

function files = download_zip(repo, url, opt)
   files = unzip(url, outputdir(opt.base, repo));
end

%--------------------------------------------------------------------------

function files = download_files(repo, url, opt)
   odir = outputdir(opt.base, repo);

   [ok, msg, id] = ensure_dir_exists(odir);

   if ~ ok
      error(id, ['Unable to ensure existence of writable ', ...
                 'output directory ''%s'': %s'], odir, msg);
   end

   if ~ iscellstr(url), url = { url }; end

   % Don't wait after last request.
   wait = [ repmat(opt.pause, [1, numel(url) - 1]), 0 ];

   % Resilience for user input.
   wait((~ isfinite(wait)) | (wait < 0)) = default_wait();

   files = {};

   wsave = mrstWebSave();

   for i = 1 : numel(url)
      files = wget(files, odir, url{i}, wsave);

      pause(wait(i));
   end
end

%--------------------------------------------------------------------------

function files = move_all_to_dest(files, repo, opt)
   nchar = 8;

   comp  = split(repo);
   rev   = opt.revision(1 : min(numel(opt.revision), nchar));
   stage = fullfile(opt.base, comp{:});
   match = dir(fullfile(stage, [comp{end}, '-', rev, '*']));

   if (numel(match) > 1) && (numel(opt.revision) > nchar)
      % Uncommon case.  Hash collision in first 'nchar' characters.  Retry
      % with full user-specified SHA-1 string.
      match = dir(fullfile(stage, [comp{end}, '-', opt.revision, '*']));
   end
   if mrstPlatform('octave')
       % We need the Octave module, otherwise different behavior for dir on
       % that platform will lead to issues here.
       require octave;
   end
   if numel(match) == 1
      % Typical case.  SHA-1 unique in first 'nchar' characters.
      odir = fullfile(stage, match.name);
   else
      error('DirMatch:Failure', ...
           ['Failed to identify GitHub repository files ', ...
            'corresponding to commit ID ''%s'' in stage/base ', ...
            'directory ''%s'''], opt.revision, stage);
   end

   assert (~isempty(odir), 'Internal logic error');

   elems = dir(odir);
   elems = { elems.name };
   skip  = strcmp(elems, '.') | strcmp(elems, '..');

   % Equivalent to FULLFILE(odir, elems(~skip)) in recent MATLAB versions.
   fileset = cellfun(@(e) fullfile(odir, e), elems(~skip), ...
                     'UniformOutput', false);

   if isempty(fileset)
      files = {};
   end

   moved_elems = move_fileset_to_dest(fileset, opt);

   assert (isempty(moved_elems) || ...
           numel(moved_elems) == numel(fileset), 'Internal Logic Error');

   if numel(moved_elems) == numel(fileset)
      % Success.  Report final fileset to caller.
      escape = @(dname) ...
         regexptranslate('escape', fileparts(fullfile(dname, '.')));

      files = files(~ cellfun(@is_directory, files));
      files = regexprep(files, escape(odir), escape(opt.dest));
   end
end

%--------------------------------------------------------------------------

function files = move_fileset_to_dest(files, opt)
   dest = opt.dest;
   [ok, msg, id] = ensure_dir_exists(dest);

   if ok
      [ok, msg] = cellfun(@(f) movefile(f, dest), files, ...
                          'UniformOutput', false);

      ok = [ ok{:} ];

      if ~ all(ok)
         report_move_failure(files(~ok), msg(~ok));
      end
   else
      warning(id, ['Failed to ensure existence of ', ...
                   'destination directory %s: %s'], dest, msg);
   end

   if ~ all(ok), files = {}; end
end

%--------------------------------------------------------------------------

function report_move_failure(fileset, msg)
   args  = [ reshape(fileset, 1, []) ; ...
             reshape(msg    , 1, []) ];

   nchar = max(cellfun('prodofsize', args(1, :)));
   args  = [ repmat({ nchar }, [ 1, size(args, 2) ]) ; args ];

   msg = sprintf('  * %*s: %s\n', args{:});

   pl = '';  if size(msg, 2) > 1, pl = 's'; end

   warning('FileMove:Failure', ...
           'Failed to move file%s to destination\n%s', pl, msg);
end

%--------------------------------------------------------------------------

function url = object_url(repo, opt)
   url = ['https://github.com/', strtrim(repo), '/'];

   if opt.all
      url = all_zip(url, opt.revision);
   else
      url = fileset(url, opt);
   end
end

%--------------------------------------------------------------------------

function url = all_zip(url, rev)
   url = [url, 'archive/', rev, '.zip'];
end

%--------------------------------------------------------------------------

function url = fileset(url, opt)
   if isempty(opt.file)
      url = '';
      return
   end

   url = strcat([url, 'blob/', opt.revision, '/'], opt.file);

   if numel(url) == 1
      url = url{1};
   end
end

%--------------------------------------------------------------------------

function [ok, msg, id] = ensure_dir_exists(d)
   [ok, msg, id] = create_if_not_exists(d);

   if ~ ok, return, end

   [ok, attr, id] = fileattrib(d);

   if ~ ok
      msg = sprintf(['Failed to obtain meta-data about ', ...
                     'directory (%s)'], attr);

   elseif ~ attr.UserWrite
      id  = 'Directory:NotWritable';
      msg = 'Directory not user writable';
   end
end

%--------------------------------------------------------------------------

function [ok, msg, id] = create_if_not_exists(odir)
   [ok, msg, id] = deal(true, '', '');

   if ~ is_directory(odir)
      [ok, msg, id] = mkdir(odir);

      if ~ ok
         msg = sprintf('Failed to create directory (%s)', msg);
      end
   end
end

%--------------------------------------------------------------------------

function files = wget(files, odir, url, wsave)
   assert (ischar(url), 'Internal Error');

   cmp = split(url);

   f = wsave(fullfile(odir, cmp{end}), url, 'raw', true);

   if ~ isempty(f)
      files = [ files, { f } ];
   end
end

%--------------------------------------------------------------------------

function odir = outputdir(basedir, repo)
   comp = split(repo);
   odir = fullfile(basedir, comp{:});
end

%--------------------------------------------------------------------------

function comp = split(s)
   comp = regexp(s, '/', 'split');
end

%--------------------------------------------------------------------------

function tf = is_directory(elem)
   if exist('isfolder', 'builtin')
      tf = isfolder(elem);
   else
      tf = isdir(elem);
   end
end
