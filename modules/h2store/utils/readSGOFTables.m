function SGOF = readSGOFTables(fid,kw, ntab,ncol)
headerLine = fgetl(fid);  % '-- SG    KRG    KROG   PCOG'
headers = strsplit(strtrim(strrep(headerLine, '--', '')));
while ~feof(fid)
    if strcmpi(strtrim(fgetl(fid)), kw)
        break;
    end
end
data = readRelPermTable(fid, kw, ntab, ncol); % 1 table, 4 columns
fclose(fid);
% Convert to struct with named fields
SGOF = cell(ntab,1);
for j =1: ntab
    for i = 1:numel(headers)
    SGOF{j}(:,i) = data{j}(:, i);  % lowercase field names
    end
end