%% converted the data processed by Homer2 into NIRS_SPM's form
%Attention:the second marker of each block should reserved in MES data when converted
%the Homer2 data into NIRS_SPM data
%written by fakun in 2019.07.12
clc;clear;

%pathn = 'D:\matlab workspace\超扫描\prepration\近红外 - 副本';
pathn = 'D:\matlab workspace\fnirs_XXX\fnirsdata3_nirs2mat'
cd(pathn);
filename_all=dir(fullfile(pathn,'*.nirs'))


for sub =1

        
        filename = filename_all(1).name
        
        load(filename,'-mat');
        
        nirs_data.oxyData = squeeze(procResult.dc(:,1,:));
        nirs_data.dxyData = squeeze(procResult.dc(:,2,:));
        nirs_data.tHbData = squeeze(procResult.dc(:,3,:));
        nirs_data.nch = size(procResult.dc,3);% channel
        nirs_data.fs = 1/0.0983;%采样率
        
        %get information of vector_onset for each condition
        nirs_data.vector_onset = zeros(size(s,1),1);
        for cond = 1:size(s,2)
            cond_vector = squeeze(s(:,cond));
            time_stimulus = find(cond_vector ==1);%get the time of marker for each condition
            for mark = 1:size(time_stimulus)
                nirs_data.vector_onset(time_stimulus(mark)) = cond;%replace the marktime with condition number
            end
        end
        converted_data_name = [filename(1:end-5),'.mat'];
        save(converted_data_name,'nirs_data');
        disp([converted_data_name,' has completed!']);
end
% movefile sub*.mat .\Homer2_wavelet_NIRS_SPM;
disp('completed all !');

