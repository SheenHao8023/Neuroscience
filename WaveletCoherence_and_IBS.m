folderPath1 = 'C:\Users\ASUS\Desktop\SCZ_tACS\data\A\ROI\EO_HC'; 
folderPath2 = 'C:\Users\ASUS\Desktop\SCZ_tACS\data\B\ROI\EO_HC'; 
fileList1 = dir(fullfile(folderPath1, 'S*.mat'));
fileList2 = dir(fullfile(folderPath2, 'S*.mat'));
coherence_matrix = zeros(12, 12);
for i = 1:length(fileList1)
    filePath1 = fullfile(folderPath1, fileList1(i).name);
    [~, name, ~] = fileparts(filePath1);
    data1=load(filePath1); 
    filePath2 = fullfile(folderPath2, fileList2(i).name);
    data2=load(filePath2); 
        for m = 1:12
            for n = 1:12
                data1_col = data1.nirsdata.oxyData(:, m);
                data2_col = data2.nirsdata.oxyData(:, n);
                min_length = min(length(data1_col), length(data2_col));
                data1_col = data1_col(1:min_length);
                data2_col = data2_col(1:min_length);
                if isempty(data1_col) || isempty(data2_col) || ~all(isfinite(data1_col)) || ~all(isfinite(data2_col))
                    coherence_matrix(m, n) = NaN;
                else
                    [wcoh,wcs,f] = wcoherence(data1_col, data2_col);
                    coherence_value = wcoh (find(0.01<=f(:,1) & f(:,1)<=0.04) , :);
                    coherence_matrix(m,n) = mean(coherence_value(:));
                end
            end
        end
    save(fullfile('C:\Users\ASUS\Desktop\SCZ_tACS\data\IBS\EO_HC', [name, '.mat']), 'coherence_matrix');
end
