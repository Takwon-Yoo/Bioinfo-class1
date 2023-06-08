import os
import sys
import time
import directory as DIR
import numpy as np
from scipy.stats import sem
import matplotlib
import matplotlib.pyplot as plt


STAR_RF = [146, 160, 182, 185, 235, 213, 230, 248, 85, 85, 90, 84, 59, 59, 52, 53]
STAR_TR = [2940, 2869, 2951, 3146, 2989, 2801, 3200, 3080, 2414, 2087, 2594, 3058, 2102, 2729, 2406, 2489, 1497,
           2443, 1818, 1746, 1560, 1437, 1448, 1217, 989, 1273, 1508, 1074, 748, 939, 648, 840]
BWA_RF = [12, 7, 13, 12, 15, 15, 16, 15, 6, 6, 6, 5, 3, 4, 3, 3]
BWA_TR = [50, 43, 45, 48, 48, 42, 49, 49, 41, 35, 44, 49, 36, 40, 38, 39, 24, 39, 30, 28, 25, 25, 25, 21,
          23, 27, 27, 22, 23, 32, 19, 26]

font = {'size': 16}
matplotlib.rc('font', **font)


def time_bar():
    star = [STAR_RF, STAR_TR]
    bwa = [BWA_RF, BWA_TR]
    star_mean = []
    star_sem = []
    bwa_mean = []
    bwa_sem = []

    for i in star:
        star_mean.append(np.mean(i))
        star_sem.append(sem(i))
    for i in bwa:
        bwa_mean.append(np.mean(i))
        bwa_sem.append(sem(i))
    labels = ['RPF-seq', 'mRNA-seq']
    x = np.arange(len(labels))
    width = 0.35

    plotdir = f'{DIR.INITIAL_FILEPATH}alignment_time.png'
    fig, ax = plt.subplots(figsize=(6, 6))
    ax.bar(x - width / 2, star_mean, width, label='STAR', yerr=star_sem)
    ax.bar(x + width / 2, bwa_mean, width, label='BWA', yerr=bwa_sem)

    ax.set_ylabel('Alignment time(sec)')
    ax.set_xticks(x)
    ax.set_xticklabels(labels)
    ax.legend()

    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)

    plt.tight_layout()
    plt.savefig(plotdir)
    plt.close()
    return None


time_bar()