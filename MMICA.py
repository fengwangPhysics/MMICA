import numpy as np
from scipy.optimize import minimize
from scipy.stats import gaussian_kde
from readCELL import readCELL
from misfunc import *
import sys

LOG2=np.log(2)

def MI(v, A, Nc, Ns, pc, pd):
    d = np.dot(v, A)
    kdeC = gaussian_kde(d[:Nc], bw_method='silverman')
    kdeD = gaussian_kde(d[Nc:], bw_method='silverman')
    pxC=kdeC(d); pxD=kdeD(d)
    px = pxC*pc + pxD*pd
    MI=(np.sum(np.log(pxC[:Nc]/px[:Nc]))+np.sum(np.log(pxD[Nc:]/px[Nc:])))/Ns/LOG2
    return -MI

def MIC(A, Nc, method='', verbose=True):
    COV=False
    if method=='COV' or method=='MICOV':
        COV=True

    No, Ns = A.shape
    Nd = Ns - Nc
    pc = float(Nc)/Ns; pd = 1.-pc

    if COV:
        C = np.cov(A)

    vlist=[]
    Entropy=-pc*np.log(pc)-pd*np.log(pd)
    # MIC
    for i in xrange(3):
        cons = [{'type': 'eq', 'fun': lambda x: np.dot(x, x)-1.}]
        if i>0:
            if COV:
                cons += [{'type': 'eq', 'fun': lambda x: np.dot(vlist, np.dot(C, x))}]
            else: 
                cons += [{'type': 'eq', 'fun': lambda x: np.dot(vlist, x)}]
        MInow=0.;
        for kk in xrange(10):
            v=np.random.rand(No)
            if method=='' or method=='COV':
                res = minimize(MI, x0=v, args=(A, Nc, Ns, pc, pd), constraints=cons, method='SLSQP', options={'disp':verbose})
            elif method=='MI' or method=='MICOV':
                res = minimize(MIall, x0=v, args=(vlist, A, Nc, Ns, pc, pd), constraints=cons, method='SLSQP', options={'disp':verbose})
            if MInow>res.fun:
                MInow=res.fun
                resNow=res
        res=resNow
        if verbose:
            print '# %d coef & MI:'%i, 
            #showinline(coef)
            print -res.fun, -MI(res.x, A, Nc, Ns, pc, pd)
        if i==0:
            vlist=res.x
            score=[-res.fun/Entropy]
        else:
            vlist=np.vstack([vlist, res.x])
            score.append(-res.fun/Entropy)
    return np.array(vlist), score

def MIall(v, vlist, A, Nc, Ns, pc, pd):
    if len(vlist)==0:
        return MI(v, A, Nc, Ns, pc, pd)
    vall=np.vstack([vlist, v])
    #print vall.shape
    d = np.dot(vall, A)
    kdeC = gaussian_kde(d[:, :Nc], bw_method='silverman')
    kdeD = gaussian_kde(d[:, Nc:], bw_method='silverman')
    pxC=kdeC(d); pxD=kdeD(d)
    px = pxC*pc + pxD*pd
    MI_=(np.sum(np.log(pxC[:Nc]/px[:Nc]))+np.sum(np.log(pxD[Nc:]/px[Nc:])))/Ns/LOG2
    return -MI_

if __name__=='__main__':
    otulist=['g__Roseburia', 'g__Actinobacillus', 'g__Aggregatibacter', 'g__[Eubacterium]', 'g__Eikenella', 'g__Lachnospira', 'g__Fusobacterium', 'g__Holdemania', 'g__Haemophilus', 'g__Coprococcus', 'g__Blautia', 'g__Campylobacter', 'g__Faecalibacterium', 'g__Ruminococcus', 'g__Neisseria', 'g__Turicibacter', 'g__Dorea', 'g__Bifidobacterium', 'g__[Ruminococcus]', 'g__Veillonella', 'g__Erwinia', 'g__Parabacteroides', 'g__Dialister', 'g__Anaerostipes', 'g__Clostridium', 'g__Bacteroides', 'g__Rothia', 'g__Porphyromonas', 'g__Phascolarctobacterium', 'g__Staphylococcus', 'g__Sutterella']
    A, Nc=readCELL(level='genus', otulist=otulist)
    if len(sys.argv)==1:
        method=""
    else:
        method=sys.argv[1]
    vlist, score=MIC(A, Nc, method=method, verbose=True)
    #vlist, score=MIC(A, Nc, method='', verbose=True)
    #vlist, score=MIC(A, Nc, method='COV', verbose=True)
    #vlist, score=MIC(A, Nc, method='MI', verbose=True)
    #vlist, score=MIC(A, Nc, method='MICOV', verbose=True)
    print ""
    d=np.dot(vlist, A)
    for i in xrange(3):
        print 'v0:',
        showarray(vlist[i,:])
        print 'd0:',
        showarray(d[i, :])
    print score
    print "Nc=", Nc
