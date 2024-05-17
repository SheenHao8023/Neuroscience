folderPath = 'C:\Users\ASUS\Desktop\SCZ_tACS\ROI\FC_SZ\cutting\HB_SCZ\'; 
files = dir(fullfile(folderPath, '*.mat')); % ʹ��dir������ȡ�ļ���������.mat�ļ�
GrangerCausality = [];
for k = 1:length(files)
    filePath = fullfile(files(k).folder, files(k).name);
    fprintf('Loading file: %s\n', filePath);
    load(filePath);
    X = nirsdata.oxyData';
    nvar = size(X, 1); %������ݾ���X�Ĵ�С����SDͨ����
    X = cca_detrend(X);
    X = cca_rm_temporalmean(X,1); %%ȥ���ƺ͹�һ��
    [bic,aic] = cca_find_model_order(X,2,30); %%Ѱ�����ŵ�ģ�ͽ���
    ret = cca_granger_regress(X,bic,1);  %����ʱ����������������BIC��׼���Ƶ�����ģ�ͽ���,����1��ʾ����F����  
    PVAL=0.05;
    [PR,q] = cca_findsignificance(ret,PVAL,1);
    GC = ret.gc;
    GC2 = GC.*PR;  %GC2��F����ĸ������������ֵ
    GCvalue = reshape(GC2, 1, 144); % ʹ��reshape����������ת��Ϊ1��144�е��������������򣨼�����������ԭ�����һ�е����ݣ�
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




