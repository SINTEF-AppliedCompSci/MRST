function [amgclpath, boostpath] = getAMGCLDependencyPaths(varargin)
    global AMGCLPATH
    global BOOSTPATH
    opt = struct('prompt', true, 'amgcl_rev', 'a551614040f0a7b793b41a4a63386675ca61d8da');
    opt = merge_options(opt, varargin{:});
    amgclpath = AMGCLPATH;
    boostpath = BOOSTPATH;
    dep_path = getDependencyFolder();
    if isempty(amgclpath)
        amgclpath = fullfile(dep_path, 'amgcl');
    end
    if isempty(boostpath)
        boostpath = fullfile(dep_path, 'boost-1_65_1_subset');
    end
    amgcl_missing = ~valid_global_path(amgclpath);
    boost_missing = ~valid_global_path(boostpath);

    if amgcl_missing
        s = sprintf('Did not find AMGCL repository path in default location "%s" or AMGCLPATH global variable. Would you like to download the files (approximately 1 MB download)?', amgclpath);
        do_download = queryDownload(s);
        if do_download
            if ~isdir(dep_path) %#ok
                mkdir(dep_path);
            end
            repo = 'ddemidov/amgcl';
            % Latest tested
            fprintf('Downloading AMGCL...')
            githubDownload(repo, 'All', true, 'Base', mrstOutputDirectory, ...
                                 'Dest', dep_path, 'Revision', opt.amgcl_rev);
            [ok, msg] = movefile(fullfile(dep_path, ['amgcl-', opt.amgcl_rev]), amgclpath);
            fprintf('Ok!');
            if ~ok
                error(msg);
            end
        end
    end
    
    if boost_missing
        s = sprintf('Did not find boost repository path in default location "%s" or BOOSTPATH global variable. Would you like to download the the requisite boost subset (approximately 1.8 MB download, may take about one minute to unzip)?', boostpath);
        do_download = queryDownload(s);
        boost_url = 'https://www.sintef.no/contentassets/124f261f170947a6bc51dd76aea66129/boost-1_65_1_subset.zip';
        if do_download
            if ~isdir(dep_path) %#ok
                mkdir(dep_path);
            end
            fprintf('Downloading and extracting BOOST. This may take a minute...')
            unzip(boost_url, dep_path);
            fprintf('Ok!');
        end
    end
end


function tf = valid_global_path(p)
   tf = ~isempty(p) && ischar(p) && isdir(p);
end

function status = queryDownload(msg)
    isDesktop = usejava('desktop');
    if isDesktop
        title = 'Missing dependency';
        choice = questdlg(msg, title,'Yes','No', 'Yes');
    else
        prompt = [msg, ' y/n [y]: '];
        choice = input(prompt,'s');
    end
    status = strcmpi(choice, 'y') || strcmpi(choice, 'yes');
end

function fldr = getDependencyFolder()
    fldr = fullfile(mrstPath('linearsolvers'), 'amgcl', 'dependencies');
end