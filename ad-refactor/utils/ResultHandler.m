classdef ResultHandler < handle
    properties
        writeToDisk
        storeInMemory
        
        dataDirectory
        dataFolder
        dataPrefix
        saveflags
        
        data
        verbose
    end
    
    methods
        
        function handler = ResultHandler(varargin)
            handler.writeToDisk = true;
            handler.storeInMemory = false;

            handler.dataDirectory = fullfile(mrstPath('query', 'ad-fi'));
            handler.dataPrefix = 'state';
            handler.dataFolder = 'cache';
            handler.saveflags = '';
            
            handler.verbose = mrstVerbose();
            
            handler = merge_options(handler, varargin{:});
            
            if ~(handler.writeToDisk || handler.storeInMemory)
                warning('ResultHandler:Noop',...
                    ['Configured without any storage targets!'...
                    ' All data will be discarded!']);
            end
            
            % Storage stuff
            handler.data = {};
            if handler.writeToDisk
                p = handler.getDataPath();
                if ~exist(p, 'dir') == 7
                    ok = mkdir(p);
                    if ~ok
                        error(['Unable to create output directory! ',...
                        'Ensure that you have write permissions to ''',...
                        handler.dataDirectory, ...
                        '''', ...
                        ])
                    end
                end
            end
            
        end
        
        function varargout = subsref(handler, s)
            switch s(1).type
                case '.'
                    % Methods and so on can be dispatched to matlab
                    [varargout{1:nargout}] = builtin('subsref', handler, s);
                case {'()', '{}'}
                    if handler.storeInMemory
                        [varargout{1:nargout}] = builtin('subsref', handler.data, s);
                        return
                    end
                    
                    if handler.writeToDisk
                        
                    end
            end
            
        end
        
        function handler = subsasgn(handler, s, v)
            switch s.type
                case '.'
                    handler = builtin('subsasgn', handler, s, v);
                case {'()', '{}'}
                    if iscell(s.subs)
                        assert(numel(s.subs) == 1);
                        s.subs = s.subs{1};
                    end
                    
                    if ~iscell(v)
                        if numel(v) == 1
                            [tmp{1:numel(s.subs)}] = deal(v);
                            v = tmp;
                        else
                            v = arrayfun(@(x) x, v, 'UniformOutput', false);
                        end
                    end
                    
                    if numel(v) ~= numel(s.subs)
                        error('ResultHandler:dimMismatch',...
                            'Left hand side and right hand side must have equal dimensions')
                    end
                    
                    if handler.writeToDisk
                        for i = 1:numel(s.subs)
                            handler.writeToFile(v{i}, s.subs(i));
                        end
                    end
                    
                    if handler.storeInMemory
                        for i = 1:numel(s.subs)
                            handler.data{s.subs(i)} = v{i};
                        end
                    end
                otherwise
                
            end
        end
        
        function handler = writeToFile(handler, data, id)
            p = handler.getDataPath(id);
            save(p, 'data', handler.saveflags);
            dispif(handler.verbose, 'Writing data to %s\n', p);
        end
        
        function p = getDataPath(handler, i)
            p = fullfile(handler.dataDirectory, handler.dataFolder);
            if nargin > 1
                p = fullfile(p, [handler.dataPrefix, num2str(i), '.mat']);
            end
        end
        
        function resetData(handler)
            handler.data = {};
            if handler.writeToDisk
                disp('Do something here')
            end
        end
    end
end