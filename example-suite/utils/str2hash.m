function hash = str2hash(str)
    try
        % Calculate hash value
        md   = java.security.MessageDigest.getInstance('SHA-256');
        hash = sprintf('%2.2x', typecast(md.digest(uint8(str)), 'uint8')');
    catch
        % ... fallback using string
        hash = str;
    end
end