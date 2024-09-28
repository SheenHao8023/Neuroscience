folderPath1 = 'C:\Users\ASUS\Desktop\SCZ_tACS\data\IBS\EO_SZ';
fileList1 = dir(fullfile(folderPath1, 'S*.mat'));
outputMatrix = [];
for i = 1:length(fileList1)
    fileName = fullfile(folderPath1, fileList1(i).name);
    load(fileName);
    row_coherence = coherence_matrix(:)';
    outputMatrix = [outputMatrix; row_coherence];
end
folderPath2 = 'C:\Users\ASUS\Desktop\SCZ_tACS\data\IBS\EO_HC';
fileList2 = dir(fullfile(folderPath2, 'S*.mat'));
for j = 1:length(fileList2)
    fileName = fullfile(folderPath2, fileList2(j).name);
    load(fileName);
    row_coherence = coherence_matrix(:)';
    outputMatrix = [outputMatrix; row_coherence];
end
selected_output = outputMatrix(:, [2:12 15:24 28:36 41:48 54:60 67:72 80:84 93:96 106:108 119:120 132]);
writematrix(selected_output, 'C:\Users\ASUS\Desktop\SCZ_tACS\data\IBS\EO.xlsx');
