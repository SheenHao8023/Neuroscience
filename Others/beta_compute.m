clc;
clear all;
cd ('D:\matlab workspace\fnirs_XXX\script')
load('condition.mat')

%% 提取数据名称
% pathn = 'D:\matlab workspace\fnirs_XXX\fnirsdata1_nirs2mat'
%pathn = 'D:\matlab workspace\fnirs_XXX\fnirsdata2_nirs2mat'
 pathn = 'D:\matlab workspace\fnirs_XXX\fnirsdata3_nirs2mat'
% pathn = 'D:\matlab workspace\fnirs_XXX\fnirsdata4_nirs2mat'
cd(pathn);
% load('mark1_all.mat')
%load('mark2_all.mat')
 load('mark3_all.mat')
% load('mark4_all.mat')


%% 确定onset

tic
for i=[1:27]%[31:33 37:45]%[49:66]%i=[1:15 19:30 34:36 46:48]

    data_name=filename_all(i).name
    load(data_name)
    hb = 'hbo';
    HPF = 'wavelet';%高通滤波方式
    LPF = 'hrf';%低通滤波方式
    method_cor = 0; % precoloring,进行
    dir_save =[data_name(1:end-4),'HbO\'];
    mkdir(dir_save)
    flag_window = 0;
    hrf_type = 2; % hrf + time +xxx
    units = 0; % scans
    
    %% 确认mark 和 duration
%     load('mark_onset.mat')
%     load('reac_end.mat')
    names{1}='1';%条件命名
    names{2}='2'
%     names{3}='3'
%     names{4}='4'
%     names{5}='5'
%     names{6}='6'
%     names{7}='7'
%     names{8}='8'
%% condtion1
%%condition
% onsets{1}=mark1_all{i,1}(3:7,1)
% onsets{2}=mark1_all{i,3}(3:7,1)
% durations{1}=mark1_all{i,2}(3:7,1)-mark1_all{i,1}(3:7,1)
% durations{2}=mark1_all{i,4}(3:7,1)-mark1_all{i,3}(3:7,1)
%%control
% onsets{1}=mark1_all{i,1}(1:2,1)
% onsets{2}=mark1_all{i,3}(1:2,1)
% durations{1}=mark1_all{i,2}(1:2,1)-mark1_all{i,1}(1:2,1)
% durations{2}=mark1_all{i,4}(1:2,1)-mark1_all{i,3}(1:2,1)
%% condition2
%condition
onsets{1}=mark3_all{i,1}(3:7,1)
onsets{2}=mark3_all{i,3}(3:7,1)
durations{1}=mark3_all{i,2}(3:7,1)-mark3_all{i,1}(3:7,1)
durations{2}=mark3_all{i,4}(3:7,1)-mark3_all{i,3}(3:7,1)
%control
% onsets{1}=mark4_all{i,1}(1:2,1)
% onsets{2}=mark4_all{i,3}(1:2,1)
% durations{1}=mark4_all{i,2}(1:2,1)-mark4_all{i,1}(1:2,1)
% durations{2}=mark4_all{i,4}(1:2,1)-mark4_all{i,3}(1:2,1)
%% condition3
%condition
%onsets{1}=mark3_all{i,1}(3:7,1)
% onsets{2}=mark3_all{i,3}(3:7,1)
% durations{1}=mark3_all{i,2}(3:7,1)-mark2_all{i,1}(3:7,1)
% durations{2}=mark3_all{i,4}(3:7,1)-mark2_all{i,3}(3:7,1)
%control
% onsets{1}=mark3_all{i,1}(3:7,1)
% onsets{2}=mark3_all{i,3}(3:7,1)
% durations{1}=mark3_all{i,2}(3:7,1)-mark2_all{i,1}(3:7,1)
% durations{2}=mark3_all{i,4}(3:7,1)-mark2_all{i,3}(3:7,1)
%% condition4
% onsets{1}=mark4_all{i,1}(1:5,1)
% onsets{2}=mark4_all{i,3}(1:5,1)
% durations{1}=mark4_all{i,2}(1:5,1)-mark4_all{i,1}(1:5,1)
% durations{2}=mark4_all{i,4}(1:5,1)-mark4_all{i,3}(1:5,1)
%control
% onsets{1}=mark4_all{i,1}(6:7,1)
% onsets{2}=mark4_all{i,3}(6:7,1)
% durations{1}=mark4_all{i,2}(6:7,1)-mark4_all{i,1}(6:7,1)
% durations{2}=mark4_all{i,4}(6:7,1)-mark4_all{i,3}(6:7,1)

    [SPM_nirs] = specification_batch(data_name, hb, HPF, LPF, method_cor, dir_save, flag_window, hrf_type, units, names, onsets, durations);
    
    %% estimation
    fname_SPM = [data_name(1:end-4),'HbO','\SPM_indiv_HbO.mat'];
    fname_nirs =[data_name];
    SPM_nirs = estimation_batch(fname_SPM, fname_nirs);%估计beta值的函数
    disp('Get beta values...')
    load('con_betalist.mat')
    betaList(i,1:44)=SPM_nirs.nirs.beta(1,:);%提取β值
    betaList(i,45:88)=SPM_nirs.nirs.beta(4,:)
%     betaList(i,38:74)=SPM_nirs.nirs.beta(4,:)
%     betaList(i,38:74)=SPM_nirs.nirs.beta(4,:)
%     betaList(i,1:37)=SPM_nirs.nirs.beta(1,:);%提取β值
%     betaList(i,38:74)=SPM_nirs.nirs.beta(4,:)
%     betaList(i,75:111)=SPM_nirs.nirs.beta(7,:);%提取β值
%     betaList(i,112:148)=SPM_nirs.nirs.beta(10,:)
%     betaList(i,149:185)=SPM_nirs.nirs.beta(13,:);%提取β值
%     betaList(i,186:222)=SPM_nirs.nirs.beta(16,:)
%     betaList(i,223:259)=SPM_nirs.nirs.beta(19,:)
%     betaList(i,260:296)=SPM_nirs.nirs.beta(22,:)

%     savename=['condition4_64_79','.mat']
%        savename=['control4_64_79','.mat'];
       
%            savename=['condition3_80_86','.mat']
%        savename=['control3_80_86','.mat'];

            savename=['condition2_47_55','.mat']
%        savename=['control2_47_55','.mat'];
      
    %save(savename,'betaList')
    save(savename,'betaList','-append')
    clear betaList
    
end
toc

 %% 提取beta
pathn = 'D:\matlab workspace\fnirs_XXX\fnirsdata3_nirs2mat'
cd(pathn);
% load('mark1_all.mat')
% load('mark2_all.mat')
 load('mark3_all.mat')
%load('mark4_all.mat')

for i=[1:27]%[31:33 37:45]%[49:66]%[1:15 19:30 34:36 46:48]
    data_name=filename_all(i).name
    load(data_name)
    fname_SPM = [data_name(1:end-4),'HbO']
    cd(fname_SPM)
    load('SPM_indiv_HbO.mat')
    betaList(i,1:44)=SPM_nirs.nirs.beta(1,:);%提取β值
    betaList(i,45:88)=SPM_nirs.nirs.beta(4,:)
    cd ..\
%     savename=['control1_1_20','.mat']
    %savename=['condition1_1_20','.mat']
%     savename=['control4_56_63','.mat']
     savename=['condition2_47_55','.mat']
%     savename=['control2_47_55','.mat']
%     savename=['condition4_56_63','.mat']
%     savename=['control4_64_79','.mat']
%     savename=['condition4_64_79','.mat']
%     savename=['control3_80_86','.mat']
%     savename=['condition3_80_86','.mat']
%         savename=['control2_74_78','.mat']
%     savename=['condition2_74_78','.mat']
    save(savename,'betaList')
    %save(savename,'betaList','-append')

end