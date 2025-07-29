from basicrta.util import simulate_hn
from basicrta.gibbs import Gibbs
import numpy as np
from scipy.optimize import linear_sum_assignment as lsa

#def test_gibbs():
#    wts = np.array([0.90, 0.09, 0.01])
#    rts = [5, 0.05, 0.001]
#    x = simulate_hn(1e3, wts, rts)
#    g = Gibbs(times=x, residue='X1', ncomp=3, niter=200)
#    g.g = 1
#    g.burnin = 100
#    
#    g._prepare()
#    #g.run()
#
#    #for i in range(len(G.results.mcrates)):
#    #    tmpsum = np.ones((3, 3), dtype=np.float64)
#    #    for ii in range(3):
#    #        for jj in range(3):
#    #            tmpsum[ii,jj] = abs(G.results.mcrates[i][ii]-rts[jj])
#
#    #    # Hungarian algorithm for minimum cost 
#    #    sortinds = lsa(tmpsum)[1]
#
#    #    # Relabel states
#    #    G.results.mcweights[i] = G.results.mcweights[i][sortinds] 
#    #    G.results.mcrates[i] = G.results.mcrates[i][sortinds]
#
#    #tmp = np.array([np.sort(G.results.rates[:,i]) for i in range(G.results.ncomp)])
#    #tmp2 = (tmp.cumsum(axis=1).T/tmp.cumsum(axis=1).T[-1])
#    #tmp3 = tmp.T[[np.where((tmp2[:,i]>0.025)&(tmp2[:,i]<0.975))[0] for i in range(G.results.ncomp)][0]]
#    #descsort = np.median(G.results.mcrates[1000:], axis=0).argsort()[::-1]
#    #ci = np.array([[line[0],line[-1]] for line in tmp3.T])
#
#    #Bools = np.array([(rts[i]>ci[descsort][i,0])&(rts[i]<ci[descsort][i,1]) for i in descsort])
#
#    #assert Bools.all() == True
#    assert len(g.t>0)
##    assert len(g.results) > 0 

def test_simdata():
    wts = np.array([0.90, 0.09, 0.01])
    rts = [5, 0.05, 0.001]
    x = simulate_hn(1e5, wts, rts)
    assert len(np.unique(x))==len(x)

def test_get_dec():
    pass

if __name__=="__main__":
    test_simdata()
#    test_gibbs()

