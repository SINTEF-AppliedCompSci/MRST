function [ s ] = updateFaultTrans( s, fault_finx, transMult )
% Transmissibilities for all fault faces are updated according to faultType
% faultType could be sealing (trans=0) or conducting (trans=trans*transMult)

%     opt.faultType = 'sealing';
%     opt.transMult = [];
%     opt = merge_options(opt, varargin{:});
    

    %% Modify transmissibility values for fault face indexes
    
    s.T_all( fault_finx ) = s.T_all( fault_finx ).*transMult;
    



end

