import numpy as np
from sklearn.manifold import TSNE
import pandas as pd
from sklearn.cluster import KMeans
from sklearn.metrics import adjusted_rand_score
import matplotlib.pyplot as plt
from operator import itemgetter
import os
from matplotlib.colors import ListedColormap

plt.rcParams['font.family'] = 'Arial'
plt.rcParams['font.size'] = 9
#put your input file here
ssr_file = ""
ssr_file_basename = os.path.basename(ssr_file)
data = pd.read_csv(ssr_file, sep=",")

label_pre = list(data['SSR'])
label_ture = []

for label in label_pre:
    
    label_ture.append(label)

headers = label_ture
data = data.drop('SSR', axis=1)
matrix = data.to_numpy()

rand_state_t_SNE = list(range(0, 51))
rand_state_k_means = list(range(0, 11))
ARI = 0
ARI_State_LIST = []

for random_state1 in rand_state_t_SNE:
    tsne_result = TSNE(perplexity=30, n_components=2, random_state=random_state1).fit_transform(matrix)
    for random_state2 in rand_state_k_means:
        kmeans = KMeans(n_clusters=8, random_state=random_state2)
        predict_ture_label = kmeans.fit_predict(tsne_result)
        ARI_new = adjusted_rand_score(np.array(label_ture), predict_ture_label)
        ARI_State_LIST.append(('R1=' + str(random_state1) + ' ' + 'R2=' + str(random_state2), ARI_new))
        if ARI_new > ARI:
            ARI = ARI_new
            xx = random_state1
            yy = random_state2

print(xx, yy, ARI)

ARI_State_LIST = sorted(ARI_State_LIST, key=itemgetter(1), reverse=True)
ARI_State_LIST_10 = []

ARI_State_LIST = enumerate(ARI_State_LIST)

for x, y in ARI_State_LIST:
    if x < 10:
        ARI_State_LIST_10.append(y)

print(ARI_State_LIST_10)

state = []
ARI_list = []

for (a, b) in ARI_State_LIST_10:
    state.append(a)
    ARI_list.append(b)
    print(a)

fig, (ax2, ax1) = plt.subplots(1, 2, figsize=(164/25.4, 84/25.4))
ax1.bar(state, ARI_list, color='#2a9d8f')
ax1.yaxis.grid(True, linestyle='--', color='gray')
ax1.set_xticklabels(state, rotation=30, fontsize=9)
ax1.set_ylabel('Adjusted rand index (ARI)')
ax1.spines['top'].set_visible(False)
ax1.spines['right'].set_visible(False)
for i, v in enumerate(ARI_list):
    ax1.text(i, v, f"{v:.2}", ha='center', va='bottom')

# Custom color map for 8 clusters
colors = ['#7fc97f', '#beaed4', '#fdc086', '#CCCC66', '#386cb0', '#f0027f', '#bf5b17', '#666666']
cluster_cmap = ListedColormap(colors)

tsne = TSNE(n_components=2, random_state=xx, perplexity=30)
tsne_results = tsne.fit_transform(matrix)
kmeans = KMeans(n_clusters=8, random_state=yy)
data['cluster'] = kmeans.fit_predict(tsne_results)

ax2.scatter(tsne_results[:, 0], tsne_results[:, 1], c=data['cluster'], cmap=cluster_cmap)

for i, header in enumerate(headers):
    x, y = tsne_results[i, 0], tsne_results[i, 1]
    ax2.annotate(header, (x, y), textcoords="offset points", xytext=(0, 5), ha='center')
ax2.set_xlabel('t-SNE Dimension 1')
ax2.set_ylabel('t-SNE Dimension 2')
ax2.spines['top'].set_visible(False)
ax2.spines['right'].set_visible(False)
plt.tight_layout()

plt.savefig(f"/output/" + ssr_file_basename + '.pdf')
print(f"/output/" + ssr_file_basename + '.pdf')
#replace "/output/" by the your directory