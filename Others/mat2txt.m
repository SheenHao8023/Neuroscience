mat_file_names = dir(fullfile('C:\Users\ASUS\Desktop\SCZ_tACS\ROI\NBS\*.mat'));
for i = 1:length(mat_file_names)
mat_file_name = fullfile('C:\Users\ASUS\Desktop\SCZ_tACS\ROI\NBS', mat_file_names(i).name);
load(mat_file_name); 
data = fillmissing(data, 'constant', nanmean(data(:))) ;
data (logical(eye(size(data)))) = 1;
[filepath,name,ext] = fileparts(mat_file_name);
writematrix(data, string(name)+'.txt' ,'Delimiter','space');
end