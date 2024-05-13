import mne
import os
import re
import numpy as np
import gdown
import h5py
import pywt
import pandas as pd
import matplotlib.pyplot as plt
import scipy.io
import seaborn as sns
from scipy.stats import ttest_1samp, ttest_rel
from mne.stats import fdr_correction, f_mway_rm
from neurora.stuff import clusterbased_permutation_1d_1samp_1sided, \
 permutation_test, \
 clusterbased_permutation_2d_1samp_2sided, \
 clusterbased_permutation_2d_2sided
from mne.time_frequency import tfr_array_morlet

# 查看epoch数据，分段为[-800ms,7000ms],epoch数不等，挑选了正确反应的trial？或者经过visual inspection删除了一些？
# raw = mne.io.read_epochs_eeglab('C:/Users/ASUS/Desktop/Graduation_Project/data/eegdata/anodal/Offline/S01-A-LX.set')
# print(raw)
# print(raw.info)

# 数据命名信息提取
pattern = re.compile(r'\.set$') # 定义正则表达式
folder_path = 'C:/Users/ASUS/Desktop/Graduation_Project/data/eegdata/anodal/Baseline/'
Anodal_Baseline = [os.path.join(folder_path, filename) for filename in os.listdir(folder_path) if pattern.search(filename)]
folder_path = 'C:/Users/ASUS/Desktop/Graduation_Project/data/eegdata/anodal/Online/'
Anodal_Online = [os.path.join(folder_path, filename) for filename in os.listdir(folder_path) if pattern.search(filename)]
folder_path = 'C:/Users/ASUS/Desktop/Graduation_Project/data/eegdata/anodal/Offline/'
Anodal_Offline = [os.path.join(folder_path, filename) for filename in os.listdir(folder_path) if pattern.search(filename)]
ABnames = []
for filename in os.listdir('C:/Users/ASUS/Desktop/Graduation_Project/data/eegdata/anodal/'):
    if filename.endswith('.h5'):
        base_name = os.path.splitext(filename)[0]
        ABnames.append(base_name)

folder_path = 'C:/Users/ASUS/Desktop/Graduation_Project/data/eegdata/sham/Baseline/'
Sham_Baseline = [os.path.join(folder_path, filename) for filename in os.listdir(folder_path) if pattern.search(filename)]
folder_path = 'C:/Users/ASUS/Desktop/Graduation_Project/data/eegdata/sham/Online/'
Sham_Online = [os.path.join(folder_path, filename) for filename in os.listdir(folder_path) if pattern.search(filename)]
folder_path = 'C:/Users/ASUS/Desktop/Graduation_Project/data/eegdata/sham/Offline/'
Sham_Offline = [os.path.join(folder_path, filename) for filename in os.listdir(folder_path) if pattern.search(filename)]
SBnames = []
for filename in os.listdir('C:/Users/ASUS/Desktop/Graduation_Project/data/eegdata/sham/'):
    if filename.endswith('.h5'):
        base_name = os.path.splitext(filename)[0]
        SBnames.append(base_name)


# 提取并存储.h5文件
for i in range(11):
    subdatafile = mne.io.read_epochs_eeglab(Anodal_Baseline[i])
    AB = subdatafile.get_data() # subdata:NumPy, shape:<72epochs,10channels,3900times
    AB = AB[:, np.delete(np.arange(AB.shape[1]), [8, 9]), :] # delete VEOG&HEOG, shape:64epochs,8channels,3900times
    # print("数据矩阵形状:", subdata.shape)
    subdatafile = mne.io.read_epochs_eeglab(Anodal_Online[i])
    AON = subdatafile.get_data() # shape:<72*3epochs,10channels,3900times
    AON = AON[:, np.delete(np.arange(AON.shape[1]), [8, 9]), :] 
    subdatafile = mne.io.read_epochs_eeglab(Anodal_Offline[i])
    AOFF = subdatafile.get_data() # shape:<72epochs,10channels,3900times
    AOFF = AOFF[:, np.delete(np.arange(AOFF.shape[1]), [8, 9]), :] 
    # 将每个被试三种条件下的脑电数据（矩阵形式）分Keys存在⼀个.h5⽂件
    f = h5py.File('C:/Users/ASUS/Desktop/Graduation_Project/data/eegdata/anodal/' + ABnames[i] + '.h5', 'w')
    f.create_dataset('Anodal_Baseline', data=AB)
    f.create_dataset('Anodal_Online', data=AON)
    f.create_dataset('Anodal_Offline', data=AOFF)
    f.close()
    # f = h5py.File('C:/Users/ASUS/Desktop/Graduation_Project/data/eegdata/anodal/' + ABnames[i] + '.h5', 'r')
    # for Anodal_Baseline in f.keys():
        #  print(f[Anodal_Baseline].name)
        #  print(f[Anodal_Baseline].shape)
for i in range(11):
    subdatafile = mne.io.read_epochs_eeglab(Sham_Baseline[i])
    SB = subdatafile.get_data() 
    SB = SB[:, np.delete(np.arange(SB.shape[1]), [8, 9]), :] 
    subdatafile = mne.io.read_epochs_eeglab(Sham_Online[i])
    SON = subdatafile.get_data()
    SON = SON[:, np.delete(np.arange(SON.shape[1]), [8, 9]), :] 
    subdatafile = mne.io.read_epochs_eeglab(Sham_Offline[i])
    SOFF = subdatafile.get_data()
    SOFF = SOFF[:, np.delete(np.arange(SOFF.shape[1]), [8, 9]), :] 
    f = h5py.File('C:/Users/ASUS/Desktop/Graduation_Project/data/eegdata/sham/' + SBnames[i] + '.h5', 'w')
    f.create_dataset('Sham_Baseline', data=SB)
    f.create_dataset('Sham_Online', data=SON)
    f.create_dataset('Sham_Offline', data=SOFF)
    f.close()

# 时频分析：时频谱分析与相干性分析 TFR power & ITC
tfr_AB = np.zeros([11, 8, 17, 3900]) # n_subjects, n_channels, n_freqs，n_times
tfr_AON = np.zeros([11, 8, 17, 3900])
tfr_AOFF = np.zeros([11, 8, 17, 3900])
itc_AB = np.zeros([11, 8, 17, 3900])
itc_AON = np.zeros([11, 8, 17, 3900])
itc_AOFF = np.zeros([11, 8, 17, 3900])
freqs = np.arange(1, 35, 2)  # 设定时频分析的参数，频段选取1-35Hz，Step=2
n_cycles = freqs / 2.
for j in range(11):
    # 读取该被试的数据
    with h5py.File('C:/Users/ASUS/Desktop/Graduation_Project/data/eegdata/anodal/' + ABnames[j] + '.h5', 'r') as f:
        subdata_AB = np.array(f['Anodal_Baseline'])
        subdata_AON = np.array(f['Anodal_Online'])
        subdata_AOFF = np.array(f['Anodal_Offline'])
        f.close()
    subtfr_AB = tfr_array_morlet(subdata_AB, 500, freqs, n_cycles, output='power')
    subtfr_AON = tfr_array_morlet(subdata_AON, 500, freqs, n_cycles, output='power')
    subtfr_AOFF = tfr_array_morlet(subdata_AOFF, 500, freqs, n_cycles, output='power')
    subitc_AB = tfr_array_morlet(subdata_AB, 500, freqs, n_cycles, output='itc')
    subitc_AON = tfr_array_morlet(subdata_AON, 500, freqs, n_cycles, output='itc')
    subitc_AOFF = tfr_array_morlet(subdata_AOFF, 500, freqs, n_cycles, output='itc')
    tfr_AB[j] = np.average(subtfr_AB, axis=0)
    tfr_AON[j] = np.average(subtfr_AON, axis=0)
    tfr_AOFF[j] = np.average(subtfr_AOFF, axis=0)
    itc_AB[j] = np.average(subitc_AB, axis=0)
    itc_AON[j] = np.average(subitc_AON, axis=0)
    itc_AOFF[j] = np.average(subitc_AOFF, axis=0)
    # 基线校正，这⾥使⽤'logratio'⽅法，即除以基线均值并取log，取基线为-400到0ms
    for chl in range(8):
        for freq in range(len(freqs)):
            tfr_AB[j, chl, freq] = 10 * np.log10(tfr_AB[j, chl, freq] /
                                    np.average(tfr_AB[j, chl, freq, 200:400]))
            tfr_AON[j, chl, freq] = 10 * np.log10(tfr_AON[j, chl, freq] /
                                    np.average(tfr_AON[j, chl, freq, 200:400]))
            tfr_AOFF[j, chl, freq] = 10 * np.log10(tfr_AOFF[j, chl, freq] /
                                    np.average(tfr_AOFF[j, chl, freq, 200:400]))
            itc_AB[j, chl, freq] = 10 * np.log10(itc_AB[j, chl, freq] /
                                    np.average(itc_AB[j, chl, freq, 200:400]))  
            itc_AON[j, chl, freq] = 10 * np.log10(itc_AON[j, chl, freq] /
                                    np.average(itc_AON[j, chl, freq, 200:400]))  
            itc_AOFF[j, chl, freq] = 10 * np.log10(itc_AOFF[j, chl, freq] /
                                    np.average(itc_AOFF[j, chl, freq, 200:400])) 

# 平均通道查看TFR结果
tfr_ave_AB = np.mean(tfr_AB, axis=1)
tfr_ave_AON = np.mean(tfr_AON, axis=1)
tfr_ave_AOFF = np.mean(tfr_AOFF, axis=1)

tfr_SB = np.zeros([11, 8, 17, 3900]) 
tfr_SON = np.zeros([11, 8, 17, 3900])
tfr_SOFF = np.zeros([11, 8, 17, 3900])
itc_SB = np.zeros([11, 8, 17, 3900])
itc_SON = np.zeros([11, 8, 17, 3900])
itc_SOFF = np.zeros([11, 8, 17, 3900])
for j in range(11):
    with h5py.File('C:/Users/ASUS/Desktop/Graduation_Project/data/eegdata/sham/' + SBnames[j] + '.h5', 'r') as f:
        subdata_SB = np.array(f['Sham_Baseline'])
        subdata_SON = np.array(f['Sham_Online'])
        subdata_SOFF = np.array(f['Sham_Offline'])
        f.close()
    subtfr_SB = tfr_array_morlet(subdata_SB, 500, freqs, n_cycles, output='power')
    subtfr_SON = tfr_array_morlet(subdata_SON, 500, freqs, n_cycles, output='power')
    subtfr_SOFF = tfr_array_morlet(subdata_SOFF, 500, freqs, n_cycles, output='power')
    subitc_SB = tfr_array_morlet(subdata_SB, 500, freqs, n_cycles, output='itc')
    subitc_SON = tfr_array_morlet(subdata_SON, 500, freqs, n_cycles, output='itc')
    subitc_SOFF = tfr_array_morlet(subdata_SOFF, 500, freqs, n_cycles, output='itc')
    tfr_SB[j] = np.average(subtfr_SB, axis=0)
    tfr_SON[j] = np.average(subtfr_SON, axis=0)
    tfr_SOFF[j] = np.average(subtfr_SOFF, axis=0)
    itc_SB[j] = np.average(subitc_SB, axis=0)
    itc_SON[j] = np.average(subitc_SON, axis=0)
    itc_SOFF[j] = np.average(subitc_SOFF, axis=0)
    for chl in range(8):
        for freq in range(len(freqs)):
            tfr_SB[j, chl, freq] = 10 * np.log10(tfr_SB[j, chl, freq] /
                                    np.average(tfr_SB[j, chl, freq, 200:400]))
            tfr_SON[j, chl, freq] = 10 * np.log10(tfr_SON[j, chl, freq] /
                                    np.average(tfr_SON[j, chl, freq, 200:400]))
            tfr_SOFF[j, chl, freq] = 10 * np.log10(tfr_SOFF[j, chl, freq] /
                                    np.average(tfr_SOFF[j, chl, freq, 200:400]))
            itc_SB[j, chl, freq] = 10 * np.log10(itc_SB[j, chl, freq] /
                                    np.average(itc_SB[j, chl, freq, 200:400]))  
            itc_SON[j, chl, freq] = 10 * np.log10(itc_SON[j, chl, freq] /
                                    np.average(itc_SON[j, chl, freq, 200:400]))  
            itc_SOFF[j, chl, freq] = 10 * np.log10(itc_SOFF[j, chl, freq] /
                                    np.average(itc_SOFF[j, chl, freq, 200:400]))  

# 平均通道查看TFR结果
tfr_ave_SB = np.mean(tfr_SB, axis=1)
tfr_ave_SON = np.mean(tfr_SON, axis=1)
tfr_ave_SOFF = np.mean(tfr_SOFF, axis=1)


# 将时频结果导出为.mat文件
scipy.io.savemat('C:/Users/ASUS/Desktop/Graduation_Project/data/eegdata/TFR.mat', {
    'tfr_AB': tfr_AB,
    'tfr_SB': tfr_SB,
    'tfr_AON': tfr_AON,
    'tfr_SON': tfr_SON,
    'tfr_AOFF': tfr_AOFF,
    'tfr_SOFF': tfr_SOFF,
})
scipy.io.savemat('C:/Users/ASUS/Desktop/Graduation_Project/data/eegdata/ITC.mat', {
    'itc_AB': itc_AB,
    'itc_SB': itc_SB,
    'itc_AON': itc_AON,
    'itc_SON': itc_SON,
    'itc_AOFF': itc_AOFF,
    'itc_SOFF': itc_SOFF,
})

# TFR_Power可视化：定义⼀个绘制单条件时频分析结果的函数
# 该函数包含统计分析与可视化的功能,并且对结果进⾏了cluster-based permutation test
def plot_tfr_results(tfr, freqs, times, p=0.01, clusterp=0.05, clim=[-4, 4]):
    n_freqs = len(freqs)
    n_times = len(times)
    # 注意：cluster-based permutation test，需要运⾏较⻓时间,iter越大越费时
    # 其返回的stats_results为⼀个shape为[n_freqs, n_times]的矩阵
    # 该矩阵中不显著的点的值为0，显著⼤于0的点的值为1，显著⼩于0的点的值为-1
    stats_results = clusterbased_permutation_2d_1samp_2sided(tfr, 0,   # NeuroRA的stuff模块
    p_threshold=p,
    clusterp_threshold=clusterp,
    iter=100)
    fig, ax = plt.subplots(1, 1)  # 时频分析结果可视化
    padsats_results = np.zeros([n_freqs + 2, n_times + 2])      # 勾勒显著性区域
    padsats_results[1:n_freqs+1, 1:n_times + 1] = stats_results
    x = np.concatenate(([times[0]-1], times, [times[-1]+1]))
    y = np.concatenate(([freqs[0]-1], freqs, [freqs[-1]+1]))
    X, Y = np.meshgrid(x, y)
    ax.contour(X, Y, padsats_results, [0.5], colors="red", alpha=0.9,linewidths=2, linestyles="dashed")
    ax.contour(X, Y, padsats_results, [-0.5], colors="blue", alpha=0.9,linewidths=2, linestyles="dashed")
    im = ax.imshow(np.average(tfr, axis=0), cmap='RdBu_r', origin='lower',
    extent=[times[0], times[-1], freqs[0], freqs[-1]], clim=clim)
    ax.set_aspect('auto')
    cbar = fig.colorbar(im)
    cbar.set_label('dB', fontsize=12)
    ax.set_xlabel('Time (ms)', fontsize=16)
    ax.set_ylabel('Frequency (Hz)', fontsize=16)
    # plt.show()    # 绘制时频结果热⼒图
    # print(stats_results)

# TFR_Power可视化：定义⼀个绘制两条件时频分析结果差异的函数
def plot_tfr_diff_results(tfr1, tfr2, freqs, times, p=0.01, clusterp=0.05, clim=[-2, 2]):
    n_freqs = len(freqs)
    n_times = len(times)
    stats_results = clusterbased_permutation_2d_2sided(tfr1, tfr2,   
    p_threshold=p,
    clusterp_threshold=clusterp,
    iter=100)
    tfr_diff = tfr1 - tfr2
    fig, ax = plt.subplots(1, 1) 
    padsats_results = np.zeros([n_freqs + 2, n_times + 2])    
    padsats_results[1:n_freqs+1, 1:n_times + 1] = stats_results
    x = np.concatenate(([times[0]-1], times, [times[-1]+1]))
    y = np.concatenate(([freqs[0]-1], freqs, [freqs[-1]+1]))
    X, Y = np.meshgrid(x, y)
    ax.contour(X, Y, padsats_results, [0.5], colors="red", alpha=0.9,linewidths=2, linestyles="dashed")
    ax.contour(X, Y, padsats_results, [-0.5], colors="blue", alpha=0.9,linewidths=2, linestyles="dashed")
    im = ax.imshow(np.average(tfr_diff, axis=0), cmap='RdBu_r', origin='lower',
    extent=[times[0], times[-1], freqs[0], freqs[-1]], clim=clim)
    ax.set_aspect('auto')
    cbar = fig.colorbar(im)
    cbar.set_label('dB', fontsize=12)
    ax.set_xlabel('Time (ms)', fontsize=16)
    ax.set_ylabel('Frequency (Hz)', fontsize=16)
    # plt.show()    # 绘制时频结果热⼒图
    # print(stats_results)


# 各电极点作图，时间范围选取[-400,2100]ms
# ch_names:k in range(8)：PZ, P3, P4, POZ, O1, O2, PO7, PO8   
times = np.arange(-400, 2100, 2)
for k in range(8)      
    tfr_AB1 = tfr_AB[:, k, :, 200:1450] #以k通道，时间范围选取[-400,2100]ms进行作图
    tfr_SB1 = tfr_SB[:, k, :, 200:1450]
    tfr_AON1 = tfr_AON[:, k, :, 200:1450]
    tfr_SON1 = tfr_SON[:, k, :, 200:1450]
    tfr_AOFF1 = tfr_AOFF[:, k, :, 200:1450]
    tfr_SOFF1 = tfr_SOFF[:, k, :, 200:1450]
    FAB = plot_tfr_results(tfr_AB1, freqs, times, p=0.05, clusterp=0.05, clim=[-3, 3])
    FSB = plot_tfr_results(tfr_SB1, freqs, times, p=0.05, clusterp=0.05, clim=[-3, 3])
    FAON = plot_tfr_results(tfr_AON1, freqs, times, p=0.05, clusterp=0.05, clim=[-3, 3])
    FSON = plot_tfr_results(tfr_SON1, freqs, times, p=0.05, clusterp=0.05, clim=[-3, 3])
    FAOFF = plot_tfr_results(tfr_AOFF1, freqs, times, p=0.05, clusterp=0.05, clim=[-3, 3])
    FSOFF = plot_tfr_results(tfr_SOFF1, freqs, times, p=0.05, clusterp=0.05, clim=[-3, 3])
    FDB = plot_tfr_diff_results(tfr_AB1, tfr_SB1, freqs, times, p=0.05, clusterp=0.05, clim=[-2, 2])
    FDON = plot_tfr_diff_results(tfr_AON1, tfr_SON1, freqs, times, p=0.05, clusterp=0.05, clim=[-2, 2])
    FDOFF = plot_tfr_diff_results(tfr_AOFF1, tfr_SOFF1, freqs, times, p=0.05, clusterp=0.05, clim=[-2, 2])
plt.show()

# 平均通道作图，时间范围选取[-400,2100]ms  
times = np.arange(-400, 2100, 2)   
tfr_AB1 = tfr_ave_AB[:, :, 200:1450] #时间范围选取[-400,2100]ms进行作图
tfr_SB1 = tfr_ave_SB[:, :, 200:1450]
tfr_AON1 = tfr_ave_AON[:, :, 200:1450]
tfr_SON1 = tfr_ave_SON[:, :, 200:1450]
tfr_AOFF1 = tfr_ave_AOFF[:, :, 200:1450]
tfr_SOFF1 = tfr_ave_SOFF[:, :, 200:1450]
# FAB = plot_tfr_results(tfr_AB1, freqs, times, p=0.05, clusterp=0.05, clim=[-3, 3])
# FSB = plot_tfr_results(tfr_SB1, freqs, times, p=0.05, clusterp=0.05, clim=[-3, 3])
# FAON = plot_tfr_results(tfr_AON1, freqs, times, p=0.05, clusterp=0.05, clim=[-3, 3])
# FSON = plot_tfr_results(tfr_SON1, freqs, times, p=0.05, clusterp=0.05, clim=[-3, 3])
# FAOFF = plot_tfr_results(tfr_AOFF1, freqs, times, p=0.05, clusterp=0.05, clim=[-3, 3])
# FSOFF = plot_tfr_results(tfr_SOFF1, freqs, times, p=0.05, clusterp=0.05, clim=[-3, 3])
FDB = plot_tfr_diff_results(tfr_AB1, tfr_SB1, freqs, times, p=0.05, clusterp=0.05, clim=[-3, 3])
FDON = plot_tfr_diff_results(tfr_AON1, tfr_SON1, freqs, times, p=0.05, clusterp=0.05, clim=[-3, 3])
FDOFF = plot_tfr_diff_results(tfr_AOFF1, tfr_SOFF1, freqs, times, p=0.05, clusterp=0.05, clim=[-3, 3])
plt.show()


