function data = ReadTextFile(filePath, SkipNoLines)

fid = fopen(filePath, 'r');
for k=1:SkipNoLines
    tline = fgets(fid);
end
formatSpec = '%f';
sizeA = [SkipNoLines-1 Inf];
data = fscanf(fid,formatSpec,sizeA);
data = data';
fclose(fid);
