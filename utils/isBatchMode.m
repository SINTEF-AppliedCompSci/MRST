function is_batch = isBatchMode()
% Check if MRST is running in batch mode
%
% SYNOPSIS:
%   is_batch = isBatchMode()
%
% NOTE:
%   This function simply checks if the global MRST_BATCH exists, and is
%   true if it exists.

    global MRST_BATCH
    if isempty(MRST_BATCH)
        MRST_BATCH = false;
    end
    assert(islogical(MRST_BATCH), 'MRST_BATCH global must be a true or false');
    is_batch = MRST_BATCH;
end
