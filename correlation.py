from ossaudiodev import SOUND_MIXER_SYNTH
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from pandas.plotting import parallel_coordinates, scatter_matrix
from math import sqrt

data = pd.read_csv('./Leaf/leaf.csv', \
                    names=["Class","Specimen Number","Eccentricity","Aspect Ratio","Elongation", "Solidity", \
                    "Stochastic Convexity", "Isoperimetric Factor", "Maximal Indentation Depth", "Lobedness", \
                    "Average Intensity", "Average Contrast", "Smoothness", "Third moment", "Uniformity", "Entropy"])

names=["Class","Eccentricity","Aspect Ratio","Elongation", "Solidity", \
        "Stochastic Convexity", "Isoperimetric Factor", "Maximal Indentation Depth", "Lobedness", \
        "Average Intensity", "Average Contrast", "Smoothness", "Third moment", "Uniformity", "Entropy"]
names_no_class =["Eccentricity","Aspect Ratio","Elongation", "Solidity", \
        "Stochastic Convexity", "Isoperimetric Factor", "Maximal Indentation Depth", "Lobedness", \
        "Average Intensity", "Average Contrast", "Smoothness", "Third moment", "Uniformity", "Entropy"]

correlation = [[0 for _ in range(len(names))] for _ in range(len(names))]

for i in range(len(names)):
    for j in range(i, len(names)):
        if i == j:
            continue

        avg_i = data[names[i]].mean()
        avg_j = data[names[j]].mean()

        sum_x_y = 0
        sum_sqr_x = 0
        sum_sqr_y = 0

        for k in range(len(data[names[i]])):
            x_i = data[names[i]][k]
            x_j = data[names[j]][k]
            sum_x_y += (x_i - avg_i)*(x_j - avg_j)
            sum_sqr_x += (x_i - avg_i)**2
            sum_sqr_y += (x_j - avg_j)**2

        correlation[i][j] = sum_x_y/(sqrt(sum_sqr_x)*sqrt(sum_sqr_y))
        correlation[j][i] = correlation[i][j]

for i in range(len(correlation)):
    for j in range(len(correlation)):
        if abs(correlation[i][j]) < 0.1:
            correlation[i][j] = 0

a = np.copy(correlation)
d = np.zeros((len(a), len(a)))
sum_corr = np.sum(a, axis = 0)

for i in range(len(d)):
    d[i, i] = sum_corr[i]

l = d-a
eig, vec = np.linalg.eig(l)
#print(eig)

eig_min, eig_idx = 10, 0

for i in range(len(eig)):
    if eig[i] < eig_min and eig[i] > 0.0001:
        eig_min = eig[i]
        eig_idx = i
#print(eig_min)

vect = [[vec[eig_idx][i], [i]] for i in range(len(a))]
vect.sort()

vect_idx = [vect[j][1] for j in range(len(vect))]
#print(vect_idx)

for _ in range(14):
    mini = abs(vect[0][0] - vect[1][0])
    idx = 0
    for i in range(len(vect)-1):
        diff = abs(vect[i][0] - vect[i+1][0])
        if diff < mini and diff > 0:
            mini = diff
            idx = i
    merge = vect[idx][1]+vect[idx+1][1]

    #print(mini, merge)

    vect[idx][1] = merge

    vect[idx][0] = 0
    for i in merge:
        vect[idx][0] += vec[eig_idx][i]
    vect[idx][0] /= len(merge)

    del vect[idx+1]
    #print()
    #print(vect)

groups = [[2, 3, 4, 5, 6, 9], [1], [7], [8], [10, 12], [11, 14], [13]]
group_names = []
for g in groups :
    for i in range(1,len(g)):
        group_names.append(names[g[i]])

normalised = data.copy()

for name in names_no_class:
    maximum = max(normalised[name])
    for i in range(len(data[name])):
        normalised[name][i] /= maximum

for group in groups :
    for i in range(len(normalised[names[0]])):
        total = 0
        for k in group :
            total += normalised[names[k]][i]
        normalised[group[0]] = total/len(group)

for name in group_names:
    del normalised[name]

choice = list(map(int,input("\nEnter the numbers : ").strip().split()))
chosen = normalised["Class"].isin(choice)

from random import randint
colors = []
for i in range(len(choice)):
    colors.append('#%06X' % randint(0, 0xFFFFFF))

plot = normalised.loc[chosen]
plot.rename(columns = {'Aspect Ratio':'C/AR/El/S/SC/IF/AI', 'Eccentricity':'E', 'Maximal Indentation Depth':'MID', 'Lobedness':'L', 'Average Contrast':'AC/TM', 'Smoothness': 'S/En', 'Uniformity':'U'}, inplace = True)

parallel_coordinates(plot, 'Class', cols=['C/AR/El/S/SC/IF/AI', 'E', 'MID', 'L', 'AC/TM', 'S/En', 'U'], color=colors)
plt.show()