classdef VFPTable
    properties
        refDepth
        flowType
        gasRatioType
        waterRatioType
    end
    
    properties (Access = private)
        % Table for BHP
        bhp_table
        % Interpolation operator for current table
        interpolator
        % Flow lookup (always present)
        flow
        % THP lookup (always present)
        thp
        
        % Additional parameter support
        parameterValues
        parameterNames
    end
    
    methods
        function table = VFPTable(varargin)
            if nargin == 1
                table = table.setTableFromDeck(varargin{1});
            else
                table = table.updateTable(varargin{:});
            end
        end
        
        function table = updateInterpolants(table)
            x = {table.flow, table.thp, table.parameterValues{:}};
            Y = table.bhp_table;
            table.interpolator = getMultiDimInterpolator(x, Y);
        end
        
        function table = updateTable(table, flow, thp, pval, pname, T, type, varargin)
            struct('refDepth', 0, ...
                   'gasRatioType', '', ...
                   'waterRatioType', '');

            table.bhp_table = T;
            table.thp = thp;
            table.flow = flow;
            table.flowType = type;
            table.parameterValues = pval;
            table.parameterNames = pname;
            
            table = merge_options(table, varargin{:});
            table = table.updateInterpolants();
        end
        
        function bhp = evaluateBHP(table, flow, thp, varargin)
            bhp = table.interpolator(flow, thp, varargin{:}); 
        end
        
        function table = setTableFromDeck(table, d)
            if isfield(d, 'Q')
                % Producer - VFPPROD
                pval = {d.WFR, d.GFR};
                pname = {'WFR', 'GFR'};
                T = d.Q;
                extra = {'waterRatioType', lower(d.WFRID),...
                         'gasRatioType',   lower(d.GFRID)};
                if numel(d.ALQ) > 1
                    if isempty(d.ALQID)
                        warning('No ALQ field specified -- but table was provided! Selecting first table.')
                        T = squeeze(T(:, :, :, :, 1));
                    else
                        extra = [extra, d.ALQID];
                    end
                end
            else
                % Injector VFPINJ
                pval = {};
                pname = {};
                T = d.BHP;
                extra = {};
            end
            table = table.updateTable(d.FLO, d.THP, pval, pname, T, lower(d.FLOID),...
                'refDepth', d.depth, extra{:});
        end
    end
end