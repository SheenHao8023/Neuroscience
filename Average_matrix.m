D = dir('C:\Users\ASUS\Desktop\SCZ_tACS\ROI analysis\FC_SZ\cutting\HA_SCZ\CORR\*.mat');
HA = zeros(12,12);
for i = 1:length(D)
load(['C:\Users\ASUS\Desktop\SCZ_tACS\ROI analysis\FC_SZ\cutting\HA_SCZ\CORR\' D(i).name]);
data( find (isnan(data)==1)) = 0.5;
HA = HA + data;
end
HA = HA / 23;