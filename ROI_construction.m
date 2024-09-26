folderPath = 'C:\Users\ASUS\Desktop\SCZ_tACS\data\A\segmentation\EO_HC'; 
filePattern = fullfile(folderPath, 'S*.mat');
fileList = dir(filePattern);
for i = 1:length(fileList)
    filePath = fullfile(folderPath, fileList(i).name);
    load(filePath); 
    zscored_oxyData = zscore(nirsdata.oxyData);
        newdata = zeros(size(zscored_oxyData, 1), 12); % 初始化为0，共有12列
        newdata(:, 1) = mean(zscored_oxyData(:, [1 2 4 8]), 2); % 第一列
        newdata(:, 2) = mean(zscored_oxyData(:, [3 6 7]), 2);   % 第二列
        newdata(:, 3) = mean(zscored_oxyData(:, [11 13 14 15]), 2); % 第三列
        newdata(:, 4) = mean(zscored_oxyData(:, [12 16]), 2);     % 第四列
        newdata(:, 5) = mean(zscored_oxyData(:, [19 20]), 2);     % 第五列
        newdata(:, 6) = mean(zscored_oxyData(:, [17 21 22 23]), 2); % 第六列
        newdata(:, 7) = mean(zscored_oxyData(:, [26 30 31 35]), 2); % 第七列
        newdata(:, 8) = mean(zscored_oxyData(:, [27 28]), 2);       % 第八列
        newdata(:, 9) = mean(zscored_oxyData(:, [36 45]), 2);       % 第九列
        newdata(:, 10) = mean(zscored_oxyData(:, [32 33 34 40]), 2); % 第十列
        newdata(:, 11) = mean(zscored_oxyData(:, [37 38 46]), 2);    % 第十一列
        newdata(:, 12) = mean(zscored_oxyData(:, [42 44 47 48]), 2); % 第十二列
    nirsdata.oxyData = newdata;
    nirsdata.nch=12;
    nirsdata.exception_channel=[zeros(1,12)];
    save(filePath, 'nirsdata');
end
