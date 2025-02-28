folderPath = 'C:\Users\ASUS\Desktop\SCZ_tES\data\ROI\FC_SZ\cutting\HB_SCZ\'; 
files = dir(fullfile(folderPath, '*.mat')); % 使用dir函数获取文件夹中所有.mat文件
GrangerCausality = [];
for k = 1:length(files)
    filePath = fullfile(files(k).folder, files(k).name);
    fprintf('Loading file: %s\n', filePath);
    load(filePath);
    X = nirsdata.oxyData';
    nvar = size(X, 1); %获得数据矩阵X的大小，即SD通道数
    X = cca_detrend(X);
    X = cca_rm_temporalmean(X,1); %%去趋势和归一化
    [bic,aic] = cca_find_model_order(X,2,30); %%寻找最优的模型阶数
    ret = cca_granger_regress(X,bic,1);  %计算时域格兰杰因果，采用BIC标准估计的最优模型阶数,参数1表示进行F检验  
    PVAL=0.05;
    [PR,q] = cca_findsignificance(ret,PVAL,1);
    GC = ret.gc;
    GC2 = GC.*PR;  %GC2是F检验的格兰杰因果连接值
    GCvalue = reshape(GC2, 1, 144); % 使用reshape函数将矩阵转换为1行144列的向量，按列排序（即向量先排列原矩阵第一列的内容）
    GrangerCausality(k,:) = GCvalue;
    %figure(6?)
    % for i=1:12
    %     for j=1:12
    %     subplot(12,12,(i-1)*12+j)
    %     stem(GC2(i,j))
    %     end
    % end
end
writematrix(GrangerCausality, 'C:/Users/ASUS/Desktop/HB_SZ.xlsx');




