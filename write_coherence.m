folderPath = 'C:\Users\ASUS\Desktop\SCZ_tACS\ROI\FC_SZ\cutting\EO_SCZ\CORR';
filePattern = fullfile(folderPath, 'EO*.mat');
fileList = dir(filePattern);
outputMatrix = [];

for i = 1:length(fileList)
    fileName = fullfile(folderPath, fileList(i).name);
    data = load(fileName);
    if isfield(data, 'data')
        rowVector = data.data(:)';
        outputMatrix = [outputMatrix; rowVector];
    end
end
selected_output = outputMatrix(:, [2:12 15:24 28:36 41:48 54:60 67:72 80:84 93:96 106:108 119:120 132]);

