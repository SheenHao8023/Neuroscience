
%将所有被试的beta值进行组分析，计算激活
cd('D:\matlab workspace\fnirs_XXX\结果')
clear all
close all
clc
cd
load('mixed_anova_cond.mat')

Alpha_Level=0.05

%
for k=1:3%condition
    for j=1:88
        if k==1
            cbeta_channel=mix_anova_cond(1:69,j+2);
        elseif k==2
            cbeta_channel=mix_anova_cond(1:69,88+j+2)
        elseif k==3
            cbeta_channel=mix_anova_cond(1:69,88*2+j+2)
        end
        
        [h,p,ci,stats]=ttest(cbeta_channel,0,'Alpha',Alpha_Level);
        % Test the null hypothesis that the sample data comes from a
        % population with mean equal to zero
        t=stats.tstat;
        %      mu=0; x=cbeta_channel; xbar = mean(x); s = std(x); t=
        %      (xbar-mu)/(s/sqrt(nSub));
        
        T_Values(j,k)=t;%通道*条件的p值
        P_Value(j,k)=p;
    end
end

% ************ FDR ********************************

% y = fdr0(p, q) to calculate whether a pvalue survive FDR corrected q p:
% an array of p values. (e.g. p values for each channel) q: desired FDR
% threshold (typically 0.05 or 0.01) y: an array of the same size with p
% with only two possible values. 0 means this position (channel) does not
% survive the threshold, 1 mean it survives
% http://www.alivelearn.net/?p=1914
y = fdr0(P_Value(:,:),0.05);

%save('3_con.mat','T_Values','P_Value','y')
save('overall_con_2nd_all.mat','T_Values','P_Value','y')

 [fdr1_x,fdr1_y]=find(y==1)
[p1_x,p1_y]=find(P_Value<0.05)
 p1=[p1_x p1_y]
 [t1_x,t1_y]=find(T_Values>0)
 t1=[t1_x t1_y]
 save('overall_con_2nd_all.mat')

% ********************** 2D T-map *********************

% plotTopoMap(T_Values(:,1)', '3x5', [-2 2]);
%  Copy_of_topomap_plot_3x7(T_Values(18:34,4)',[-3,3])






