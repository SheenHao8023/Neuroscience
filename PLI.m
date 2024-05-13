%  Anodal条件
h5Files = dir(fullfile('C:/Users/ASUS/Desktop/Graduation_Project/data/eegdata/anodal/', '*.h5'));   % MATLAB是时间-通道-epoch顺序
fullFilePathsWithFullPath = cell(size(h5Files));  % 创建一个cell数组来存储完整路径
for i = 1:length(h5Files)
    fullFilePathsWithFullPath{i} = fullfile('C:/Users/ASUS/Desktop/Graduation_Project/data/eegdata/anodal/', h5Files(i).name);
end
start_sample = 600;
end_sample = 1200;  % 定义采样点范围，时间对应400-1600ms
result_matrix = zeros(1, 64);  % 建立空矩阵用于PLI指标统计
for sub = 1:11
    eeg = h5read(fullFilePathsWithFullPath{sub}, '/Anodal_Baseline');
    subset_eeg = eeg(start_sample:end_sample, :, :);  % 提取指定采样点范围内的数据
    nevents = size(subset_eeg, 3); % epoch的数量
    total_pli = zeros(8, 8, nevents); % 初始化用于存储每个epoch的8×8 PLI矩阵
    avg_pli = zeros(8, 8, 11); % 初始化用于存储每个被试的8×8 PLI矩阵
    for epoch_idx = 1:nevents % 遍历每个epoch
        epoch_data = subset_eeg(:, :, epoch_idx); % 提取当前epoch的数据
        for ch1 = 1:8
            for ch2 = ch1+1:8 % 避免重复计算对称项（PLI(ch1, ch2) = PLI(ch2, ch1)）
                thetax=angle(hilbert(epoch_data(:, ch1)));
                thetay=angle(hilbert(epoch_data(:, ch2)));
                thetadiff=thetax-thetay;
                thetadiffmodify=wrapToPi(thetadiff);
                thetadiffmodify((find(thetadiffmodify)==pi)|(find(thetadiffmodify)==-pi))=0;
                PLI=abs(mean(sign(thetadiffmodify)));    % 计算PLI
                total_pli(ch1, ch2, epoch_idx) = PLI;
                total_pli(ch2, ch1, epoch_idx) = PLI;  % 由于无向FC是对称的，所以也可以直接赋值给对角线位置
            end
        end
    end
    avg_pli = mean(total_pli, 3); % 第三维度上求平均值，平均所有的epoch得到每个被试最终的功能连接矩阵
    save(sprintf('C:/Users/ASUS/Desktop/Graduation_Project/data/eegdata/anodal/Baseline/S%02d.mat', sub) , 'avg_pli');
    result_matrix(sub,:) = avg_pli(:);
end
condtion_pli = mean(avg_pli, 3); % 第三维度上求平均值，平均所有的被试得到最终的功能连接矩阵
save('C:/Users/ASUS/Desktop/Graduation_Project/data/eegdata/AB.mat', 'condtion_pli');
writematrix(result_matrix, 'C:/Users/ASUS/Desktop/Graduation_Project/data/eegdata/AB.xlsx');

for sub = 1:11
    eeg = h5read(fullFilePathsWithFullPath{sub}, '/Anodal_Online');
    subset_eeg = eeg(start_sample:end_sample, :, :);  % 提取指定采样点范围内的数据
    nevents = size(subset_eeg, 3); % epoch的数量
    total_pli = zeros(8, 8, nevents); % 初始化用于存储每个epoch的8×8 PLI矩阵
    avg_pli = zeros(8, 8, 11); % 初始化用于存储每个被试的8×8 PLI矩阵
    for epoch_idx = 1:nevents % 遍历每个epoch
        epoch_data = subset_eeg(:, :, epoch_idx); % 提取当前epoch的数据
        for ch1 = 1:8
            for ch2 = ch1+1:8 % 避免重复计算对称项（PLI(ch1, ch2) = PLI(ch2, ch1)）
                thetax=angle(hilbert(epoch_data(:, ch1)));
                thetay=angle(hilbert(epoch_data(:, ch2)));
                thetadiff=thetax-thetay;
                thetadiffmodify=wrapToPi(thetadiff);
                thetadiffmodify((find(thetadiffmodify)==pi)|(find(thetadiffmodify)==-pi))=0;
                PLI=abs(mean(sign(thetadiffmodify)));    % 计算PLI
                total_pli(ch1, ch2, epoch_idx) = PLI;
                total_pli(ch2, ch1, epoch_idx) = PLI;  % 由于无向FC是对称的，所以也可以直接赋值给对角线位置
            end
        end
    end
    avg_pli = mean(total_pli, 3); % 第三维度上求平均值，平均所有的epoch得到每个被试最终的功能连接矩阵
    save(sprintf('C:/Users/ASUS/Desktop/Graduation_Project/data/eegdata/anodal/Online/S%02d.mat', sub) , 'avg_pli');
    result_matrix(sub,:) = avg_pli(:);
end
condtion_pli = mean(avg_pli, 3); % 第三维度上求平均值，平均所有的被试得到最终的功能连接矩阵
save('C:/Users/ASUS/Desktop/Graduation_Project/data/eegdata/AON.mat', 'condtion_pli');
writematrix(result_matrix, 'C:/Users/ASUS/Desktop/Graduation_Project/data/eegdata/AON.xlsx');

for sub = 1:11
    eeg = h5read(fullFilePathsWithFullPath{sub}, '/Anodal_Offline');
    subset_eeg = eeg(start_sample:end_sample, :, :);  % 提取指定采样点范围内的数据
    nevents = size(subset_eeg, 3); % epoch的数量
    total_pli = zeros(8, 8, nevents); % 初始化用于存储每个epoch的8×8 PLI矩阵
    avg_pli = zeros(8, 8, 11); % 初始化用于存储每个被试的8×8 PLI矩阵
    for epoch_idx = 1:nevents % 遍历每个epoch
        epoch_data = subset_eeg(:, :, epoch_idx); % 提取当前epoch的数据
        for ch1 = 1:8
            for ch2 = ch1+1:8 % 避免重复计算对称项（PLI(ch1, ch2) = PLI(ch2, ch1)）
                thetax=angle(hilbert(epoch_data(:, ch1)));
                thetay=angle(hilbert(epoch_data(:, ch2)));
                thetadiff=thetax-thetay;
                thetadiffmodify=wrapToPi(thetadiff);
                thetadiffmodify((find(thetadiffmodify)==pi)|(find(thetadiffmodify)==-pi))=0;
                PLI=abs(mean(sign(thetadiffmodify)));    % 计算PLI
                total_pli(ch1, ch2, epoch_idx) = PLI;
                total_pli(ch2, ch1, epoch_idx) = PLI;  % 由于无向FC是对称的，所以也可以直接赋值给对角线位置
            end
        end
    end
    avg_pli = mean(total_pli, 3); % 第三维度上求平均值，平均所有的epoch得到每个被试最终的功能连接矩阵
    save(sprintf('C:/Users/ASUS/Desktop/Graduation_Project/data/eegdata/anodal/Offline/S%02d.mat', sub) , 'avg_pli');
    result_matrix(sub,:) = avg_pli(:);
end
condtion_pli = mean(avg_pli, 3); % 第三维度上求平均值，平均所有的被试得到最终的功能连接矩阵
save('C:/Users/ASUS/Desktop/Graduation_Project/data/eegdata/AOFF.mat', 'condtion_pli');
writematrix(result_matrix, 'C:/Users/ASUS/Desktop/Graduation_Project/data/eegdata/AOFF.xlsx');

%  Sham条件
h5Files = dir(fullfile('C:/Users/ASUS/Desktop/Graduation_Project/data/eegdata/sham/', '*.h5'));   % MATLAB是时间-通道-epoch顺序
fullFilePathsWithFullPath = cell(size(h5Files));  % 创建一个cell数组来存储完整路径
for i = 1:length(h5Files)
    fullFilePathsWithFullPath{i} = fullfile('C:/Users/ASUS/Desktop/Graduation_Project/data/eegdata/sham/', h5Files(i).name);
end
start_sample = 600;
end_sample = 1200;  % 定义采样点范围，时间对应400-1600ms
result_matrix = zeros(1, 64);  % 建立空矩阵用于PLI指标统计
for sub = 1:11
    eeg = h5read(fullFilePathsWithFullPath{sub}, '/Sham_Baseline');
    subset_eeg = eeg(start_sample:end_sample, :, :);  % 提取指定采样点范围内的数据
    nevents = size(subset_eeg, 3); % epoch的数量
    total_pli = zeros(8, 8, nevents); % 初始化用于存储每个epoch的8×8 PLI矩阵
    avg_pli = zeros(8, 8, 11); % 初始化用于存储每个被试的8×8 PLI矩阵
    for epoch_idx = 1:nevents % 遍历每个epoch
        epoch_data = subset_eeg(:, :, epoch_idx); % 提取当前epoch的数据
        for ch1 = 1:8
            for ch2 = ch1+1:8 % 避免重复计算对称项（PLI(ch1, ch2) = PLI(ch2, ch1)）
                thetax=angle(hilbert(epoch_data(:, ch1)));
                thetay=angle(hilbert(epoch_data(:, ch2)));
                thetadiff=thetax-thetay;
                thetadiffmodify=wrapToPi(thetadiff);
                thetadiffmodify((find(thetadiffmodify)==pi)|(find(thetadiffmodify)==-pi))=0;
                PLI=abs(mean(sign(thetadiffmodify)));    % 计算PLI
                total_pli(ch1, ch2, epoch_idx) = PLI;
                total_pli(ch2, ch1, epoch_idx) = PLI;  % 由于无向FC是对称的，所以也可以直接赋值给对角线位置
            end
        end
    end
    avg_pli = mean(total_pli, 3); % 第三维度上求平均值，平均所有的epoch得到每个被试最终的功能连接矩阵
    save(sprintf('C:/Users/ASUS/Desktop/Graduation_Project/data/eegdata/sham/Baseline/S%02d.mat', sub) , 'avg_pli');
    result_matrix(sub,:) = avg_pli(:);
end
condtion_pli = mean(avg_pli, 3); % 第三维度上求平均值，平均所有的被试得到最终的功能连接矩阵
save('C:/Users/ASUS/Desktop/Graduation_Project/data/eegdata/SB.mat', 'condtion_pli');
writematrix(result_matrix, 'C:/Users/ASUS/Desktop/Graduation_Project/data/eegdata/SB.xlsx');

for sub = 1:11
    eeg = h5read(fullFilePathsWithFullPath{sub}, '/Sham_Online');
    subset_eeg = eeg(start_sample:end_sample, :, :);  % 提取指定采样点范围内的数据
    nevents = size(subset_eeg, 3); % epoch的数量
    total_pli = zeros(8, 8, nevents); % 初始化用于存储每个epoch的8×8 PLI矩阵
    avg_pli = zeros(8, 8, 11); % 初始化用于存储每个被试的8×8 PLI矩阵
    for epoch_idx = 1:nevents % 遍历每个epoch
        epoch_data = subset_eeg(:, :, epoch_idx); % 提取当前epoch的数据
        for ch1 = 1:8
            for ch2 = ch1+1:8 % 避免重复计算对称项（PLI(ch1, ch2) = PLI(ch2, ch1)）
                thetax=angle(hilbert(epoch_data(:, ch1)));
                thetay=angle(hilbert(epoch_data(:, ch2)));
                thetadiff=thetax-thetay;
                thetadiffmodify=wrapToPi(thetadiff);
                thetadiffmodify((find(thetadiffmodify)==pi)|(find(thetadiffmodify)==-pi))=0;
                PLI=abs(mean(sign(thetadiffmodify)));    % 计算PLI
                total_pli(ch1, ch2, epoch_idx) = PLI;
                total_pli(ch2, ch1, epoch_idx) = PLI;  % 由于无向FC是对称的，所以也可以直接赋值给对角线位置
            end
        end
    end
    avg_pli = mean(total_pli, 3); % 第三维度上求平均值，平均所有的epoch得到每个被试最终的功能连接矩阵
    save(sprintf('C:/Users/ASUS/Desktop/Graduation_Project/data/eegdata/sham/Online/S%02d.mat', sub) , 'avg_pli');
    result_matrix(sub,:) = avg_pli(:);
end
condtion_pli = mean(avg_pli, 3); % 第三维度上求平均值，平均所有的被试得到最终的功能连接矩阵
save('C:/Users/ASUS/Desktop/Graduation_Project/data/eegdata/SON.mat', 'condtion_pli');
writematrix(result_matrix, 'C:/Users/ASUS/Desktop/Graduation_Project/data/eegdata/SON.xlsx');

for sub = 1:11
    eeg = h5read(fullFilePathsWithFullPath{sub}, '/Sham_Offline');
    subset_eeg = eeg(start_sample:end_sample, :, :);  % 提取指定采样点范围内的数据
    nevents = size(subset_eeg, 3); % epoch的数量
    total_pli = zeros(8, 8, nevents); % 初始化用于存储每个epoch的8×8 PLI矩阵
    avg_pli = zeros(8, 8, 11); % 初始化用于存储每个被试的8×8 PLI矩阵
    for epoch_idx = 1:nevents % 遍历每个epoch
        epoch_data = subset_eeg(:, :, epoch_idx); % 提取当前epoch的数据
        for ch1 = 1:8
            for ch2 = ch1+1:8 % 避免重复计算对称项（PLI(ch1, ch2) = PLI(ch2, ch1)）
                thetax=angle(hilbert(epoch_data(:, ch1)));
                thetay=angle(hilbert(epoch_data(:, ch2)));
                thetadiff=thetax-thetay;
                thetadiffmodify=wrapToPi(thetadiff);
                thetadiffmodify((find(thetadiffmodify)==pi)|(find(thetadiffmodify)==-pi))=0;
                PLI=abs(mean(sign(thetadiffmodify)));    % 计算PLI
                total_pli(ch1, ch2, epoch_idx) = PLI;
                total_pli(ch2, ch1, epoch_idx) = PLI;  % 由于无向FC是对称的，所以也可以直接赋值给对角线位置
            end
        end
    end
    avg_pli = mean(total_pli, 3); % 第三维度上求平均值，平均所有的epoch得到每个被试最终的功能连接矩阵
    save(sprintf('C:/Users/ASUS/Desktop/Graduation_Project/data/eegdata/sham/Offline/S%02d.mat', sub) , 'avg_pli');
    result_matrix(sub,:) = avg_pli(:);
end
condtion_pli = mean(avg_pli, 3); % 第三维度上求平均值，平均所有的被试得到最终的功能连接矩阵
save('C:/Users/ASUS/Desktop/Graduation_Project/data/eegdata/SOFF.mat', 'condtion_pli');
writematrix(result_matrix, 'C:/Users/ASUS/Desktop/Graduation_Project/data/eegdata/SOFF.xlsx');

%  heatmap，Anodal/Sham两组被试分别在Baseline/Online/Offline条件下的通道间两两PLI功能连接值
figure()
SHM1=SHeatmap(condtion_pli,'Format','sq');
SHM1=SHM1.draw();
colormap(othercolor('RdBu8'))
SHM1.setText();
ax=gca;
ax.XTickLabel={'PZ', 'P3', 'P4', 'POZ', 'O1', 'O2', 'PO7', 'PO8'};
ax.YTickLabel={'PZ', 'P3', 'P4', 'POZ', 'O1', 'O2', 'PO7', 'PO8'};
ax.FontSize=18;
clim([0,0.2])
