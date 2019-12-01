import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import seaborn as sns
sns.set_style('darkgrid')

cluster_ffit=open('./cluster_fit.txt','r')
cluster_fit=np.loadtxt(cluster_ffit)
cluster_fit=np.transpose(cluster_fit)
avg_cluster_fit=np.mean(cluster_fit,0)
max_cluster_fit=np.amax(cluster_fit,0)
cluster_ffit.close()

cluster_fsig=open('./cluster_sig.txt','r')
cluster_sig=np.loadtxt(cluster_fsig)
cluster_sig=np.transpose(cluster_sig)
avg_cluster_sig=np.mean(cluster_sig,0)
cluster_fsig.close()

full_ffit=open('./full_fit.txt','r')
full_fit=np.loadtxt(full_ffit)
full_fit=np.transpose(full_fit)
avg_full_fit=np.mean(full_fit,0)
max_full_fit=np.amax(full_fit,0)
full_ffit.close()

full_fsig=open('./full_sig.txt','r')
full_sig=np.loadtxt(full_fsig)
full_sig=np.transpose(full_sig)
avg_full_sig=np.mean(full_sig,0)
full_fsig.close()




fig=plt.figure()
gs=gridspec.GridSpec(2,1)

ax_fit=fig.add_subplot(gs[0,:])
ax_fit.plot(avg_full_fit,label='full_avg')
ax_fit.plot(max_full_fit,label='full_max')
ax_fit.plot(avg_cluster_fit,label='cluster_avg')
ax_fit.plot(max_cluster_fit,label='cluster_max')
ax_fit.legend()

ax_sig=fig.add_subplot(gs[1,:])
#ax_sig.plot(avg_sig,label='avg')
ax_sig.boxplot(cluster_sig)
ax_sig.legend()
plt.show()
