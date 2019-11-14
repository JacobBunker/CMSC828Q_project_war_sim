import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import seaborn as sns
sns.set_style('darkgrid')

ffit=open('./fit.txt','r')
fit=np.loadtxt(ffit)
fit=np.transpose(fit)
avg_fit=np.mean(fit,0)
max_fit=np.amax(fit,0)

ffit.close()

fsig=open('./sig.txt','r')
sig=np.loadtxt(fsig)
sig=np.transpose(sig)
avg_sig=np.mean(sig,0)
fsig.close()



fig=plt.figure()
gs=gridspec.GridSpec(2,1)

ax_fit=fig.add_subplot(gs[0,:])
ax_fit.plot(avg_fit,label='avg')
ax_fit.plot(max_fit,label='max')
ax_fit.legend()

ax_sig=fig.add_subplot(gs[1,:])
#ax_sig.plot(avg_sig,label='avg')
ax_sig.boxplot(sig)
ax_sig.legend()
plt.show()
