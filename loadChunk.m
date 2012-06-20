function snd = loadChunk(dirPath)

files = dir(dirPath);

fInd = randi(length(files));
name = files(fInd).name;
while files(fInd).isdir | ~strcmp(name(end-3:end), '.wav')
    fInd = randi(length(files));
    name = files(fInd).name;
end

snd = wavread(strcat(dirPath, name));

fprintf(sprintf('%s loaded\n', name));