for i = 1:12
for j = 1:12
[wcoh,wcs,f] = wcoherence(nirsdata.oxyData(:,i) , nirsdata.oxyData(:,j));
A = wcoh (find(0.01<=f(:,1) & f(:,1)<=0.04) , :);
data(i,j) = mean(A(:));
end
end