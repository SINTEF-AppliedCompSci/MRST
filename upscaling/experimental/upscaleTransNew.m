function [HT_cg, T_cg, cgwells, report] = upscaleTransNew(cg, T_fine, varargin)
    warning('This interface is outdated! Use upscaleTrans directly.');
    [HT_cg, T_cg, cgwells, report] = upscaleTrans(cg, T_fine, varargin{:});
end