import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from matplotlib.font_manager import FontProperties

# Create DataFrames
Intron = pd.DataFrame({
    'Intron1': [284188, 281536, 281717, 283140, 286881, 288406, 280281, 279751, 284278, 281293, 284780, 286727, 287011, 0, 290478, 216380],
    'Intron2': [180544, 178089, 186512, 181722, 0, 175231, 180575, 179918, 180515, 184087, 183908, 174045, 174436, 0, 192020, 141566],
    'Intron3': [61272, 63250, 61825, 62751, 0, 63039, 62670, 62856, 63481, 63478, 64383, 61908, 63321, 0, 71222, 59917],
    'Intron4': [146956, 149215, 150698, 150436, 0, 152535, 143221, 142242, 144179, 155926, 156484, 157275, 160535, 0, 166883, 133878],
    'Intron5': [80673, 78327, 77920, 78442, 0, 78688, 78383, 78412, 78599, 81166, 81538, 75865, 74976, 0, 76013, 57828],
    'Intron6': [187393, 187521, 188497, 187740, 188842, 188890, 187547, 186366, 188701, 191791, 187487, 196862, 197744, 0, 195021, 130930],
    'Inron7': [216355, 214819, 215939, 217183, 208811, 208827, 217333, 216726, 219575, 215926, 223293, 216978, 215656, 0, 228957, 197561],
    'Intron8': [20351, 20336, 20303, 18443, 20108, 20104, 18788, 18992, 18705, 18955, 18363, 20472, 20854, 10279, 19920, 8437],
    'Intron9': [161976, 165829, 163344, 163343, 162114, 162242, 162549, 162311, 164042, 163718, 165700, 163261, 161522, 140494, 182642, 139192],
    'Intron10': [26588, 26654, 26270, 26204, 27353, 27861, 26698, 26259, 27709, 26867, 25775, 27458, 26696, 19233, 27565, 16133],
    'Intron11': [9876, 9336, 9572, 10355, 9504, 9491, 9014, 9036, 9008, 8964, 8942, 9222, 9900, 9849, 9243, 7708]
})

Exon = pd.DataFrame({
    'Exon1': [105, 151, 155, 242, 115, 118, 126, 142, 104, 142, 213, 142, 106, 0, 129, 115],
    'Exon2': [164, 164, 764, 264, 164, 164, 164, 164, 164, 164, 164, 164, 164, 0, 164, 164],
    'Exon3': [241, 241, 241, 241, 0, 241, 241, 241, 241, 61, 241, 241, 241, 0, 241, 242],
    'Exon4': [122, 122, 122, 122, 0, 122, 122, 122, 122, 122, 122, 122, 122, 0, 119, 122],
    'Exon5': [84, 84, 84, 84, 0, 84, 84, 84, 84, 84, 84, 84, 84, 0, 84, 84],
    'Exon6': [116, 116, 116, 116, 116, 116, 116, 116, 116, 116, 116, 116, 116, 0, 116, 116],
    'Exon7': [137, 137, 137, 137, 137, 137, 137, 137, 139, 137, 137, 137, 137, 0, 137, 137],
    'Exon8': [62, 62, 62, 62, 62, 62, 62, 62, 62, 62, 62, 62, 62, 773, 62, 62],
    'Exon9': [150, 150, 150, 150, 150, 150, 150, 150, 150, 150, 150, 150, 150, 150, 150, 148],
    'Exon10': [84, 84, 84, 84, 84, 84, 84, 84, 84, 84, 84, 84, 84, 84, 84, 84],
    'Exon11': [118, 118, 118, 118, 118, 118, 118, 118, 118, 118, 119, 117, 118, 118, 118, 118],
    'Exon12': [2795, 2773, 2787, 2703, 2709, 2766, 2774, 2773, 2890, 2831, 2753, 113, 2790, 1754, 2184, 2391]
})

# Calculate statistics
mean_introns = Intron.mean()
error_introns_std = Intron.std()
n_values_introns = len(Intron)
error_introns_se = error_introns_std / np.sqrt(n_values_introns)

mean_exons = Exon.mean()
error_exons_std = Exon.std()
n_values_exons = len(Exon)
error_exons_se = error_exons_std / np.sqrt(n_values_exons)

# Create plotting DataFrames
df_introns_plot = pd.DataFrame({
    'Gen': [f'Intron {i}' for i in range(1, len(Intron.columns) + 1)],
    'Promedio': mean_introns.values,
    'Error': error_introns_se.values
})

df_exons_plot = pd.DataFrame({
    'Gen': [f'Exon {i}' for i in range(1, len(Exon.columns) + 1)],
    'Promedio': mean_exons.values,
    'Error': error_exons_se.values
})

# Set up the plots
plt.style.use('default')
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(20, 10))
font = FontProperties(family='Times New Roman')

# Intron plot
sns.barplot(x='Gen', y='Promedio', data=df_introns_plot, ax=ax1, color='lightblue')
ax1.errorbar(x=range(len(df_introns_plot)), y=df_introns_plot['Promedio'],
             yerr=df_introns_plot['Error'], fmt='none', c='black', capsize=5)
ax1.set_title('Average Intron Size', fontproperties=font, fontweight='bold')
ax1.set_xlabel('Introns', fontproperties=font)
ax1.set_ylabel('Average Size (pb)', fontproperties=font)
ax1.tick_params(axis='x', rotation=45)
ax1.yaxis.set_major_formatter(plt.FuncFormatter(lambda x, p: format(int(x), ',')))
ax1.set_facecolor('white')
ax1.grid(False)

# Exon plot with secondary y-axis for Exon 12
ax2_twin = ax2.twinx()

# Plot Exons 1-11
sns.barplot(x='Gen', y='Promedio', data=df_exons_plot.iloc[:-1], ax=ax2, color='lightgreen')
ax2.errorbar(x=range(len(df_exons_plot)-1), y=df_exons_plot['Promedio'].iloc[:-1],
             yerr=df_exons_plot['Error'].iloc[:-1], fmt='none', c='black', capsize=5)

# Plot Exon 12 on secondary y-axis
sns.barplot(x=['Exon 12'], y=[df_exons_plot['Promedio'].iloc[-1]], ax=ax2_twin, color='lightgreen')
ax2_twin.errorbar(x=[11], y=[df_exons_plot['Promedio'].iloc[-1]],
                  yerr=[df_exons_plot['Error'].iloc[-1]], fmt='none', c='black', capsize=5)

ax2.set_title('Average Exon Size', fontproperties=font, fontweight='bold')
ax2.set_xlabel('Exons', fontproperties=font)
ax2.set_ylabel('Average Size (pb) for Exons 1-11', fontproperties=font)
ax2_twin.set_ylabel('Average Size (pb) for Exon 12', fontproperties=font)
ax2.tick_params(axis='x', rotation=45)
ax2.set_facecolor('white')
ax2.grid(False)

# Set y-axis limits
ax2.set_ylim(0, 300)
ax2_twin.set_ylim(0, 3000)

# Adjust layout and save
plt.tight_layout()
plt.savefig('Figure1.tiff', dpi=300, format='tiff', bbox_inches='tight')
plt.show()
