load('S*****.mat')
nirsdata.oxyData = fillmissing(nirsdata.oxyData, 'constant', nanmean(nirsdata.oxyData(:))) ;
nirsdata.oxyData = zscore(nirsdata.oxyData, 1);
A = (nirsdata.oxyData(:,1) + nirsdata.oxyData(:,2) + nirsdata.oxyData(:,4) + nirsdata.oxyData(:,8))/4;
B = (nirsdata.oxyData(:,3) + nirsdata.oxyData(:,6) + nirsdata.oxyData(:,7))/3;
C = (nirsdata.oxyData(:,11) + nirsdata.oxyData(:,13) + nirsdata.oxyData(:,14) + nirsdata.oxyData(:,15))/4;
D = (nirsdata.oxyData(:,12) + nirsdata.oxyData(:,16))/2;
E = (nirsdata.oxyData(:,19) + nirsdata.oxyData(:,20))/2;
F = (nirsdata.oxyData(:,17) + nirsdata.oxyData(:,21) + nirsdata.oxyData(:,22) + nirsdata.oxyData(:,23))/4;
G = (nirsdata.oxyData(:,26) + nirsdata.oxyData(:,30) + nirsdata.oxyData(:,31) + nirsdata.oxyData(:,35))/4;
H = (nirsdata.oxyData(:,27) + nirsdata.oxyData(:,28) )/2;
I = (nirsdata.oxyData(:,36) + nirsdata.oxyData(:,45))/2;
J = (nirsdata.oxyData(:,32) + nirsdata.oxyData(:,33) + nirsdata.oxyData(:,34) + nirsdata.oxyData(:,40))/4;
K = (nirsdata.oxyData(:,37) + nirsdata.oxyData(:,38) + nirsdata.oxyData(:,46))/3;
L = (nirsdata.oxyData(:,42) + nirsdata.oxyData(:,44) + nirsdata.oxyData(:,47) + nirsdata.oxyData(:,48))/4;
newdata = [A,B,C,D,E,F,G,H,I,J,K,L];
nirsdata.oxyData=newdata;
nirsdata.nch=12;
nirsdata.exception_channel=[zeros(1,12)];
save('S*****.mat', 'nirsdata')