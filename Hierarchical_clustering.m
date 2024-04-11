tiqu = @(x)(split(x,":"));
shuju = @(x)(str2double(x{end,1}));
cd([ "C:\Users\ASUS\Desktop\SCZ_tACS\behav\HC dyads\HC dyads\"]);
subfolders = dir('.\*'); 
corrmatrix = []

for i = 1:numel(subfolders)
        subfolder_path = fullfile("C:\Users\ASUS\Desktop\SCZ_tACS\behav\HC dyads\HC dyads\", subfolders(i).name);
        cd(subfolder_path);
        trials = dir(['.\2_*_*_*.mat']);  % condition = 2 (each other)
        pair=table;
    
    for n = 1:2:numel(trials)
        idx = (n-1)/2+1;
        k = findstr(trials(n).name,'_');
        a = trials(n).name(1:k(2));
        b = trials(n+1).name(1:k(2));

        % 跳过缺失配对的按键并记录
        if  ~strcmp(a,b)
            display(trials(n).name+" missing its pair, now j = "+j)
            return  % 删除不配对的试次
        end

        k = findstr(trials(n).name,'_'); % define the trial
        a = trials(n).name(1:k(2)); % condi & trial
        load(trials(n).name); load(trials(n+1).name);

        if isempty(x) || isempty(y)
            display("empty mat, now j = "+j)
            return
        end

        pair.trial(idx) = str2double(a(3:end-1));

        y{end} = [];
        y = y(cellfun(@isempty,y)==0);
        y = cellfun(@string,y);
        x{end} = [];
        x = x(cellfun(@isempty,x)==0);
        x = cellfun(@string,x);
        y = cellfun(tiqu,y,'UniformOutput',false);
        x = cellfun(tiqu,x,'UniformOutput',false);

        B.RT = cellfun(shuju,y);
        A.RT = cellfun(shuju,x);
        A.RT(isnan(A.RT))=[];B.RT(isnan(A.RT))=[];
        B.IOI = diff(B.RT);
        A.IOI = diff(A.RT);
        B.RT=cumsum(B.IOI);
        A.RT=cumsum(A.IOI);

        % 统计outliers
        if A.RT(1)<250 || B.RT(1)<250
            A.IOI(A.RT<250)=[]; B.IOI(B.RT<250)=[];
            B.RT=cumsum(B.IOI);A.RT=cumsum(A.IOI);
        end
        
        % 首：前8下没跟上，则trial直接删除
        if A.RT(1)>2000 || B.RT(1)>2000
            pair.deletion8(idx)=true;
        else
            pair.deletion8(idx)=false;
        end

        % 尾: 在这里把数据截取成一样的长度
        if numel(B.RT)<numel(A.RT) %B比A短
            A.IOI(numel(B.RT)+1:end)=[];
            A.RT(numel(B.RT)+1:end)=[];
        else
            B.IOI(numel(A.RT)+1:end)=[];%B比A长
            B.RT(numel(A.RT)+1:end)=[];
        end

        points = numel(A.RT);

        pair.Aoutlier(idx,1)  = sum(nonzeros(A.IOI<=median(A.IOI)-0.5*median(A.IOI)));
        pair.Aoutlier(idx,2) = sum(nonzeros(A.IOI>=median(A.IOI)+0.5*median(A.IOI)));
        pair.Boutlier(idx,1) = sum(nonzeros(B.IOI<=median(B.IOI)-0.5*median(B.IOI)));
        pair.Boutlier(idx,2) = sum(nonzeros(B.IOI>=median(B.IOI)+0.5*median(B.IOI)));

        %%%%%%%%%%%%%%%%%%
        %%%%%  Index  %%%%%%
        %%%%%%%%%%%%%%%%%%

        % cross-correlation
        [acorr, lags] = xcorr(A.IOI, B.IOI, 'normalized');
        corrmatrix{i, 1} = acorr(lags == 1);    % 提取A比B延迟1个样本点的互相关值
        corrmatrix{i, 2} = acorr(lags == 0);    % 提取滞后为0的互相关值
        corrmatrix{i, 3} = acorr(lags == -1);   % 提取B比A延迟1个样本点的互相关值
        % Pearson correlation and p-value
        [corrmatrix{i, 4}, corrmatrix{i, 7}] = corr(A.IOI(2:end), B.IOI(1:end-1));
        [corrmatrix{i, 5}, corrmatrix{i, 8}]  = corr(A.IOI, B.IOI);
        [corrmatrix{i, 6}, corrmatrix{i, 9}]  = corr(A.IOI(1:end-1), B.IOI(2:end));


    end
end

% 相关性数据提取
data = cell2mat (corrmatrix(3:42, 4:6));
% 估计最佳cluster个数
% evaCH = evalclusters(data,"linkage","CalinskiHarabasz","KList",1:10);
% evaDB = evalclusters(data,"linkage","DaviesBouldin","KList",1:10);
% evaSI = evalclusters(data,"linkage","silhouette","KList",1:10);
% 计算样本间的欧氏距离
dist_matrix = pdist(data, 'euclidean');
D_squared = squareform(dist_matrix); % 转换为平方形式
% 构建层次聚类链接矩阵（使用ward法）
linkage_matrix = linkage(dist_matrix, ['ward']);
% 可视化层次聚类树
F1 = dendrogram(linkage_matrix, 0);
% 结果整理
cluster_labels = [];
cluster_labels(:,1) = cluster(linkage_matrix,'MaxClust', 6);
% for i = [2:10]
%     cluster_labels(:,i-1) = cluster(linkage_matrix,'MaxClust', i);
% end
% cluster_labels(1:10, 10) = evaCH.CriterionValues';
% cluster_labels(1:10, 11) = evaDB.CriterionValues';
% cluster_labels(1:10, 12) = evaSI.CriterionValues';
% % figure index
% x = 2:10;
% y = cluster_labels(2:10, 11);
% F2=bar(x,y);
% 
% target_rows = find(cluster_labels(:,5) == 6);
% selected_data = data(target_rows, 1:3);
% size(target_rows, 1)
% meanlag = [mean(selected_data(:,1)), mean(selected_data(:,2)), mean(selected_data(:,3))];
% xlabel=2:3:8; %设置柱状图横坐标
% b = bar(xlabel, meanlag, 0.5); % 绘制柱状图函数包括三个参数：柱状图横坐标、柱状图数据、柱状图宽度
% ylim([-0.5 1]);
% set(gca,'XTickLabel',{'lagBA','lag0','lagAB'},'FontSize',12,'FontName','宋体'); %X轴标记设置
% set(get(gca, 'Xlabel'),'Fontname','Times New Roman','FontWeight','bold','Fontsize',18); % X轴标签字体，字号设置
% set(get(gca, 'Ylabel'),'Fontname','Times New Roman','FontWeight','bold','Fontsize',16); % Y轴标签字体，字号设置
% 



