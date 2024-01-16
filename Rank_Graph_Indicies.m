data = readtable('DC_RS.xlsx'); 
data = table2array(data);
[rows, cols] = size(data);
rankMatrix = zeros(rows, cols);
for r = 1:rows
    sortedRow = sort(data(r,:), 'descend');
    rankTemp = 1:cols;
    [~, index] = ismember(data(r,:), sortedRow);
    rankMatrix(r,:) = rankTemp(index);
end
writematrix(rankMatrix, 'DC_RS.xlsx');