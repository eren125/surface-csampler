# %%
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

# %%
df = pd.read_csv('results_100_surf.csv')
df['Structures'] = df['structure_file'].str.strip('out/').str.strip('.cif')
df_widom = pd.read_csv('Screening_CoReMOF_Dataset.csv')

# %%
# check nans
print(len(df.dropna()))
print(len(df))

df.dropna(inplace=True)
# df.fillna(33,inplace=True)
# %%
df_merge = pd.merge(df_widom,df, on="Structures", how="left")
df_plot = df_merge[~(df_merge['DISORDER']=='DISORDER')]
# %%
def split_deviation(x):
    i=0
    if np.isnan(x):
        return np.nan
    elif x>10:
        return "]10;40["
    while i < 4:
        if 10-2*(i+1)<x<=10-2*i:
            return "]%.1f;%.1f]"%(10-2*(i+1),10-2*i)
        i += 1
    return "[0;%.1f["%(10-2*(i))
def pairplot(x, y, x_label, y_label,left=-70,right=160,bottom=-70,top=160):
    hue = 'D_i'
    df_plot[hue+"_cat"] = df_plot[hue].apply(split_deviation)
    plt.figure(figsize=(10,6))
    hue_order=[split_deviation(11-2*i) for i in range(6)]
    palette = sns.color_palette("tab10",n_colors=len(hue_order))

    sns.scatterplot(data=df_plot.dropna(subset=[x,y,hue]),x=x,y=y,hue=hue+"_cat",hue_order=hue_order,palette=palette,alpha=0.8,s=10)
    plt.xlim(left=left,right=right)
    plt.ylim(bottom=bottom,top=top)

    plt.xlabel(x_label)
    plt.ylabel(y_label)
    plt.legend(title=r"$D_i$")
    plt.tight_layout()
    plt.savefig("plot/%s_vs_%s_overview"%(x,y),dpi=280)

# %%
x = 'H_Xe_0_widom'
x_label = r"$\Delta_{ads}H_0^{Xe}$"
y = 'enthalpy'
y_label = r"Boltzmann average of voronoi energies $E_{voro-B}^{Xe}$"
pairplot(x, y, x_label, y_label,left=-70,right=160,bottom=-70,top=160)

# %%
df_plot.dropna(subset=[x,y])[[x,y]].corr(method="pearson")
# %%
df_plot['time'].describe()
# %%
plt.xlim(right=30)
df_plot['time'].hist(bins=500)
# %%
df_plot['dr'] = abs(df_plot['enthalpy']-df_plot['H_Xe_0_widom'])/min(df_plot['enthalpy']-df_plot['H_Xe_0_widom'])

# %%
