import numpy as np
from scipy.stats import binom
G=100; N=G//2; i=np.arange(G+1); xg=i/G
def P(s):
    xs=np.clip(xg+s*xg*(1-xg),0,1); return binom.pmf(i[None,:],G,xs[:,None])
P0=P(0.0)
for B in (0.2,1.0,3.0):
    s=B/(4*N); Ps=P(s)
    for Ton in (N, 4*N):
        # mutation cohorts arising every generation during [0,Ton); selection acts while t<Ton; then neutral to absorption
        excess=0.0; pred=0.0
        d_on=np.zeros(G+1); d_on[1]=1  # we accumulate by linearity: sum over cohorts of distributions
        dist_sel=np.zeros(G+1); dist_neu=np.zeros(G+1); hetmass=0.0
        for t in range(Ton):
            dist_sel=dist_sel@Ps; dist_neu=dist_neu@P0
            dist_sel[1]+=1; dist_neu[1]+=1         # new cohort enters at 1/G (added, then propagated next step)
            hetmass+=dist_neu@(xg*(1-xg))
        for t in range(40*N): dist_sel=dist_sel@P0; dist_neu=dist_neu@P0
        ex=dist_sel@xg-dist_neu@xg
        print("B=%.1f Ton=%4d  exact excess=%.5f  first-order s*sum E x(1-x)=%.5f  ratio=%.3f  | exact-constant-B check psi: %.5f"%(B,Ton,ex,s*hetmass,ex/(s*hetmass), Ton*(1/G)*(B/(1-np.exp(-B))-1)))
