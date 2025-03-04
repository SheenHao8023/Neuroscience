%% MVGC: state-space method.
% L. Barnett and A. K. Seth, "Granger causality for state-space models", Phys. Rev. E 91(4) Rapid Communication, 2015

ntrials   = 1;     % number of trials
nobs      = size(nirsdata.oxyData, 1);   % number of observations per trial
nvars      = size(nirsdata.oxyData, 2);   % number of variables
regmode   = 'LWR';  % or 'OLS'
icregmode = 'LWR';  % or 'OLS'
morder    = 'BIC';  % or 'actual', 'AIC', numerical value)
momax     = 20;     % maximum model order for model order estimation
acmaxlags = 100;   % maximum autocovariance lags (empty for automatic calculation)
tstat     = 'F';    % statistical test for MVGC:  'F' for Granger's F-test (default) or 'chi2' for Geweke's chi2 test
alpha     = 0.05;   % significance level for significance test
mhtc      = 'FDR'; % 'NONE' 'BONFERRONI' 'SIDAK' 'HOLM' 'FDR' 'FDRD'

[AIC,BIC,moAIC,moBIC] = tsdata_to_infocrit(nirsdata.oxyData',momax,icregmode);  % Model order estimation
if     strcmpi(morder,'actual')
    morder = amo;
    fprintf('\nusing actual model order = %d\n',morder);
elseif strcmpi(morder,'AIC')
    morder = moAIC;
    fprintf('\nusing AIC best model order = %d\n',morder);
elseif strcmpi(morder,'BIC')
    morder = moBIC;
    fprintf('\nusing BIC best model order = %d\n',morder);
else
    fprintf('\nusing specified model order = %d\n',morder);
end
[A,SIG] = tsdata_to_var(nirsdata.oxyData', morder, regmode);  % VAR model estimation
[F,pval] = var_to_pwcgc(A, SIG, nirsdata.oxyData', regmode, tstat);  % Granger causality calculation: time domain 
sig = significance(pval, alpha, mhtc);
%figure(2); clf;
%sgtitlex('Pairwise-conditional Granger Causality - Time domain');
%subplot(1,3,1);
%plot_pw(F);
%title('Pairwise-conditional GC');
%subplot(1,3,2);
%plot_pw(pval);
%title(['p-values (' tstat '-test)']);
%subplot(1,3,3);
%plot_pw(sig);
%title(['Significant at \alpha = ' num2str(alpha)]);
