% Read the file as cell string line by line:
fid = fopen(FileName, 'r');
if fid < 0, error('Cannot open file: %s', FileName); end
Data = textscan(fid, '%s', 'delimiter', '\n', 'whitespace', '');
fclose(fid);
% Remove empty lines:
C    = deblank(Data{1});   % [EDITED]: deblank added
C(cellfun('isempty', C)) = [];
% Write the cell string:
fid = fopen(FileName, 'w');
if fid < 0, error('Cannot open file: %s', FileName); end
fprintf(fid, '%s\n', C{:});
fclose(fid);