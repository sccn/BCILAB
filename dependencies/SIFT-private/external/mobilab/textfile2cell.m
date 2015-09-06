function linesInCell = textfile2cell(file)
linesInCell = cell(1000,1);
count = 1;
fid = fopen(file,'r');
while ~feof(fid)
    linesInCell{count} = fgetl(fid);
    count = count+1;
end
fclose(fid);
I = cellfun(@isempty,linesInCell);
linesInCell(I) = [];
end