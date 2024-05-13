% 获取当前目录下的所有 .mat 文件
matFiles = dir(fullfile(pwd, '*.mat'));

% 遍历每一个 .mat 文件
for i = 1:length(matFiles)
    % 获取完整文件名和去掉扩展名的文件名
    matFilePath = fullfile(pwd, matFiles(i).name);
    baseFileName = matFiles(i).name(1:end-4); % 去掉 .mat 扩展名

    % 加载 .mat 文件中的 avg_pli 变量
    vars = load(matFilePath, 'avg_pli');
    data = vars.avg_pli;

    % 将 avg_pli 以8x8矩阵的形式写入到同名 .txt 文件中，使用空格作为分隔符
    txtFilePath = fullfile(pwd, [baseFileName '.txt']);
    fid = fopen(txtFilePath, 'w');
    for row = 1:size(data, 1)
        fprintf(fid, '%f ', data(row, :));
        fprintf(fid, '\n');
    end
    fclose(fid);
    disp(['Converted ' matFiles(i).name ' to ' baseFileName '.txt']);
end