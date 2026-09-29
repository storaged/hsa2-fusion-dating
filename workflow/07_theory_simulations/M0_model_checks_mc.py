import numpy as np
rng=np.random.default_rng(1)
psi=lambda B: B/(1-np.exp(-B)) if B!=0 else 1.0
# (1) Prop 1: log Phi vs <B> - <B>Var/12
for B in ([0.5]*3+[0.0]*2, [2,2,0,0,0,0], [1.0,0.2,0.0]):
    B=np.array(B,float); lp=np.log(np.mean([psi(b) for b in B])/np.mean([psi(-b) for b in B]))
    print("P1  logPhi=%.5f  approx=%.5f  <B>=%.3f"%(lp, B.mean()-B.mean()*B.var()/12, B.mean()))
# WF helpers (N diploids, 2N genes), vectorised over replicate alleles
N=200; G=2*N
def wf(p0, s_of_t, T, reps):
    x=np.full(reps,p0); het=np.zeros(reps); xsum=np.zeros(reps); tfix=np.full(reps,-1)
    for t in range(T):
        s=s_of_t(t); xs=x+s*x*(1-x)
        live=(x>0)&(x<1); het+=np.where(live,x*(1-x),0); xsum+=np.where(live,x,0)
        x=rng.binomial(G,np.clip(xs,0,1))/G
        tfix=np.where((x==1)&(tfix<0),t+1,tfix)
    return x,het,xsum,tfix
# (2) E int x(1-x) dt from p=1/2N  -> 2N p(1-p) ~ 1
x,het,xs,tf=wf(1/G,lambda t:0.0,20*N,400000)
print("P2  E int x(1-x) = %.3f (theory %.3f)"%(het.mean(), 2*N*(1/G)*(1-1/G)))
f=tf>0
print("C2  fixed: mean T_fix-T_orig = %.0f (4N=%d);  E int x dt | fix = %.0f (2N=%d);  eta=0.5 -> T_off-T_fix = %.0f"%(tf[f].mean(),4*N,xs[f].mean(),2*N,xs[f].mean()))
# (3) Lag proposition: selection on for [0,T1), off afterwards; derived-count excess per mutation
#     should equal s*E int_{on} x(1-x) (no extra delay after switch-off)
s=0.5/(4*N)   # B=0.5
T1=2*N; Tend=T1+6*N; reps=300000
def run(on_until):
    # new mutations arising uniformly in [0,T1) ; approximate by staggered starts: simulate one cohort per start time block
    tot=0
    for start in range(0,T1,N//2):
        x,_,_,_=wf(1/G,lambda t:(s if t+start<on_until else 0.0),Tend-start,reps//4)
        tot+=x.mean()
    return tot
sel=run(T1); neu=run(0)
print("P3  mean derived freq excess at end, with selection then drift: %.5f ; first-order prediction ~ s*(sum over cohorts E int_on x(1-x)) computed below"%((sel-neu)))
pred=0
for start in range(0,T1,N//2):
    _,het,_,_=wf(1/G,lambda t:0.0,T1-start,reps//4); pred+=s*het.mean()
print("P3  prediction: %.5f"%pred)
