SCZ = xlsread('C:\Users\ASUS\Desktop\SCZ_tACS\ROI\GrangerCausality\table.xlsx', 'Sheet1');

figure()
SHM1=SHeatmap(SCZ,'Format','sq');
SHM1=SHM1.draw();
colormap(jet)
SHM1.setText();
ax=gca;
ax.XTickLabel={'rA','rTPJ','rM','rDLPFC','rSFC','rFPC','lFPC','lSFC','lDLPFC','lM','lTPJ','lA'};
ax.YTickLabel={'rA','rTPJ','rM','rDLPFC','rSFC','rFPC','lFPC','lSFC','lDLPFC','lM','lTPJ','lA'};
ax.FontSize=14;
clim([0.3,0.7])
