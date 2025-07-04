import os
import gc
import pickle
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
from numpy.random import default_rng
from tqdm import tqdm
from MDAnalysis.analysis.base import Results
from basicrta.util import confidence_interval
from multiprocessing import Pool, Lock
import MDAnalysis as mda
from basicrta import istarmap

gc.enable()
mpl.rcParams['pdf.fonttype'] = 42
rng = default_rng()


class ParallelGibbs(object):
    """
    A module to take a contact map and run Gibbs samplers for each residue.

    :param contacts: Contact pickle file (`contacts-{cutoff}.pkl`).
    :type contacts: str
    :param nproc: Number of processes to use in running Gibbs samplers.
    :type nproc: int
    :param ncomp: Number of mixture components to use in the Gibbs sampler.
    :type ncomp: int
    :param niter: Number of iterations of the Gibbs sampler to perform.
    :type niter: int
    """

    def __init__(self, contacts, nproc=1, ncomp=15, niter=110000):
        self.cutoff = float(contacts.strip('.pkl').split('/')[-1].split('_')
                            [-1])
        self.niter = niter
        self.nproc = nproc
        self.ncomp = ncomp
        self.contacts = contacts

    def run(self, run_resids=None, g=100):
        """
        The :meth:`run` method executes the Gibbs samplers for all residues of
        `sel1` present in the contact map, or a list of resids can be provided.

        :param run_resids: Resid(s) for which to run a Gibbs sampler. 
        :type run_resids: int or list, optional
        """
        from basicrta.util import run_residue

        with open(self.contacts, 'r+b') as f:
            contacts = pickle.load(f)

        # Check if this is a combined contact file
        metadata = contacts.dtype.metadata
        is_combined = metadata and 'n_trajectories' in metadata and metadata['n_trajectories'] > 1
        if is_combined:
            print(f"WARNING: Using combined contact file with {metadata['n_trajectories']} trajectories.")
            print("WARNING: Kinetic clustering is not yet supported for combined contacts.")
            print("WARNING: The Gibbs sampler will pool all residence times together.")

        protids = np.unique(contacts[:, 0])
        if not run_resids:
            run_resids = protids

        if not isinstance(run_resids, (list, np.ndarray)):
            run_resids = [run_resids]

        rg = contacts.dtype.metadata['ag1'].residues
        resids = rg.resids
        reslets = np.array([mda.lib.util.convert_aa_code(name) for name in
                            rg.resnames])
        residues = np.array([f'{reslet}{resid}' for reslet, resid in
                             zip(reslets, resids)])
        times = [contacts[contacts[:, 0] == i][:, 3] for i in
                 run_resids]
        inds = np.array([np.where(resids == resid)[0][0] for resid in
                         run_resids])
        residues = residues[inds]
        input_list = [[residues[i], times[i].copy(), i % self.nproc,
                       self.ncomp, self.niter, self.cutoff, g, is_combined] for i in
                      range(len(residues))]

        del contacts, times
        gc.collect()

        with (Pool(self.nproc, initializer=tqdm.set_lock,
                   initargs=(Lock(),)) as p):
            try:
                for _ in tqdm(p.istarmap(run_residue, input_list),
                              total=len(residues), position=0,
                              desc='overall progress'):
                    pass
            except KeyboardInterrupt:
                pass


class Gibbs(object):
    r"""Gibbs sampler to estimate parameters of an exponential mixture for a set
    of data. Results are stored in :class:`gibbs.results`, which uses
    :class:`MDAnalysis.analysis.base.Results()`. If 'results=None' the gibbs
    sampler has not been executed, which requires calling :meth:`run`.

    :param times: Set of residence times to analyze
    :type times: array, optional
    :param residue: Residue name associated with the set of residence times
    :type residue: str
    :param loc: Used for progress bar in parallel applications
    :type loc: int
    :param ncomp: Number of exponential components to use in the mixture model
    :type ncomp: int
    :param niter: Number of iterations to run the Gibbs sampler
    :type niter: int
    :param cutoff: Cutoff value used in contact analysis, used to determine
                   directory to load/save results. Allows for multiple cutoffs
                   to be tested in directory containing contacts.
    :type cutoff: float

    EXAMPLE
    -------
    >>> from basicrta.gibbs import Gibbs
    >>> from basicrta.tests.datafiles import times
    >>> g = Gibbs(times=times, residue='W313', cutoff=7.0)
    >>> g.run()
    >>> g.process_gibbs()
    >>> g.estimate_tau()
    [1, 2, 3]

    To load a Gibbs sampler that has already been executed use the :meth:`load`
    method

    >>> g = Gibbs().load('results.pkl')

    The Gibbs sampler can be executed using the :meth:`run` method without
    processing the resulting data. Once the :meth:`process_gibbs` method is
    called, the :attr:`Gibbs.results.processed_results` attribute will be
    populated.
    """

    def __init__(self, times=None, residue=None, loc=0, ncomp=15, niter=110000,
                 cutoff=None, g=50, burnin=10000, gskip=2):
        self.times = times
        self.residue = residue
        self.niter = niter
        self.loc = loc
        self.ncomp = ncomp
        self.g = g
        self.gskip = gskip
        self.burnin = burnin
        self.cutoff = cutoff
        self.processed_results = Results()
        self._noise_cutoff = 0.4

        if times is not None:
            diff = (np.sort(times)[1:]-np.sort(times)[:-1])
            try:
                self.ts = diff[diff != 0][0]
            except IndexError:
                self.ts = times.min()
        else:
            self.ts = None

        self.keys = {'times', 'residue', 'loc', 'ncomp', 'niter', 'g', 'burnin',
                     'processed_results', 'ts', 'mcweights', 'mcrates', 't',
                     's', 'cutoff', 'indicator'}

    def __getitem__(self, item):
        return getattr(self, item)

    def _prepare(self):
        from basicrta.util import get_s
        self.t, self.s = get_s(self.times, self.ts)

        # initialize arrays
        self.indicator = np.zeros(((self.niter + 1) // self.g,
                                  self.times.shape[0]), dtype=np.uint8)
        self.mcweights = np.zeros(((self.niter + 1) // self.g, self.ncomp))
        self.mcrates = np.zeros(((self.niter + 1) // self.g, self.ncomp))

        # guess hyperparameters
        self.whypers = np.ones(self.ncomp) / [self.ncomp]
        self.rhypers = np.ones((self.ncomp, 2)) * [1, 3]

    def run(self):
        r"""
        Execute the Gibbs sampler and save the raw data to the instance of
        :class:`Gibbs`.
        """
        # initialize weights and rates
        self._prepare()
        if not os.path.exists(f'basicrta-{self.cutoff}/{self.residue}'):
            os.makedirs(f'basicrta-{self.cutoff}/{self.residue}')

        inrates = 0.5 * 10 ** np.arange(-self.ncomp + 2, 2, dtype=float)
        tmpw = 9 * 10 ** (-np.arange(1, self.ncomp + 1, dtype=float))
        weights, rates = tmpw / tmpw.sum(), inrates[::-1]

        # gibbs sampler
        for j in tqdm(range(1, self.niter+1),
                      desc=f'{self.residue}-K{self.ncomp}',
                      position=self.loc, leave=False):

            # compute probabilities (equation 7)
            tmp = weights*rates*np.exp(np.outer(-rates, self.times)).T
            psample = (tmp.T/tmp.sum(axis=1)).T

            # sample indicator
            z = np.argmax(rng.multinomial(1, psample), axis=1)

            # get indicator for each data point
            inds = [np.where(z == i)[0] for i in range(self.ncomp)]

            # compute total time and number of point for each component
            Ns = np.array([len(inds[i]) for i in range(self.ncomp)])
            Ts = np.array([self.times[inds[i]].sum() for i in range(self.ncomp)])

            # sample posteriors (equations 8 and 9)
            weights = rng.dirichlet(self.whypers+Ns)
            rates = rng.gamma(self.rhypers[:, 0]+Ns, 1/(self.rhypers[:, 1]+Ts))

            # save every g steps
            if j % self.g == 0:
                ind = j//self.g-1
                self.mcweights[ind], self.mcrates[ind] = weights, rates
                self.indicator[ind] = z

        self.save()

    def cluster(self, method="GaussianMixture", **kwargs):
        r"""
        Cluster the processed results using the methods available in
        :class:`sklearn.mixture`

        :param method: Mixture method to use
        :type method: str
        """
        # Check if this Gibbs result was created from combined contact data
        if hasattr(self, '_from_combined_contacts') and self._from_combined_contacts:
            print("INFO: Using combined contact data for clustering. "
                  "Trajectory source information is pooled together.")
            
        from sklearn import mixture
        from scipy import stats

        clu = getattr(mixture, method)
        burnin_ind = self.burnin // self.g
        data_len = len(self.times)
        wcutoff = 10 / data_len

        weights = self.mcweights[burnin_ind::self.gskip]
        rates = self.mcrates[burnin_ind::self.gskip]
        lens = np.array([len(row[row > wcutoff]) for row in weights])
        lmode = stats.mode(lens).mode
        train_param = lmode

        train_inds = np.where(lens == train_param)[0]
        train_weights = (weights[train_inds][weights[train_inds] > wcutoff].
                         reshape(-1, train_param))
        train_rates = (rates[train_inds][weights[train_inds] > wcutoff].
                       reshape(-1, train_param))

        inds = np.where(weights > wcutoff)
        aweights, arates = weights[inds], rates[inds]
        data = np.stack((aweights, arates), axis=1)

        tweights, trates = train_weights.flatten(), train_rates.flatten()
        train_data = np.stack((tweights, trates), axis=1)

        r = clu(**kwargs)
        r.fit(np.log(train_data))
        all_labels = r.predict(np.log(data))

        if self.indicator is not None:
            indicator = self.indicator[burnin_ind::self.gskip]
        else:
            indicator = self._sample_indicator()

        pindicator = np.zeros((self.times.shape[0], lmode))
        for j in np.unique(inds[0]):
            mapinds = all_labels[inds[0] == j]
            for i, indx in enumerate(inds[1][inds[0] == j]):
                tmpind = np.where(indicator[j] == indx)[0]
                pindicator[tmpind, mapinds[i]] += 1

        pindicator = (pindicator.T / pindicator.sum(axis=1)).T
        setattr(self.processed_results, 'indicator', pindicator)
        setattr(self.processed_results, 'labels', all_labels)

    def process_gibbs(self, show=True):
        r"""
        Process the samples collected from the Gibbs sampler.
        :meth:`process_gibbs` can be called multiple times to check the
        robustness of the results.
        """
        from basicrta.util import mixture_and_plot
        from scipy import stats

        data_len = len(self.times)
        wcutoff = 10/data_len
        burnin_ind = self.burnin//self.g
        inds = np.where(self.mcweights[burnin_ind::self.gskip] > wcutoff)
        indices = (np.arange(self.burnin, self.niter + 1, self.g*self.gskip)
                   [inds[0]] // self.g)
        weights = self.mcweights[burnin_ind::self.gskip] 
        rates = self.mcrates[burnin_ind::self.gskip]
        fweights, frates = weights[inds], rates[inds]

        lens = [len(row[row > wcutoff]) for row in
                self.mcweights[burnin_ind::self.gskip]]
        lmode = stats.mode(lens).mode

        self.cluster(n_init=117, n_components=lmode)
        labels, presorts = mixture_and_plot(self, show=show)
        setattr(self.processed_results, 'labels', labels)
        setattr(self.processed_results, 'indicator',
                self.processed_results.indicator[:, presorts])

        attrs = ["weights", "rates", "ncomp", "residue", "iteration", "niter"]
        values = [fweights, frates, lmode, self.residue, indices, self.niter]
        for attr, val in zip(attrs, values):
            setattr(self.processed_results, attr, val)

        self._estimate_params()
        self.save()

    def result_plot(self, remove_noise=False, **kwargs):
        """
        Generate the combined result plot with option to change kwargs without
        re-clustering.

        :param remove_noise: Option to remove noise clusters
        :type remove_noise: bool
        """
        from basicrta.util import mixture_and_plot
        mixture_and_plot(self, remove_noise=remove_noise, **kwargs)

    def _sample_indicator(self):
        indicator = np.zeros(((self.niter+1)//(self.g*self.gskip), 
                              self.times.shape[0]), dtype=np.uint8)
        burnin_ind = self.burnin//self.g
        for i, (w, r) in enumerate(zip(self.mcweights, self.mcrates)):
            # compute probabilities
            probs = w*r*np.exp(np.outer(-r, self.times)).T
            z = (probs.T/probs.sum(axis=1)).T

            # sample indicator
            s = np.argmax(rng.multinomial(1, z), axis=1)
            indicator[i] = s
        setattr(self, 'indicator', indicator)
        return indicator[burnin_ind::self.gskip]

    def save(self):
        """
        Save current state of the :class:`Gibbs` instance.
        """
        savedir = f'basicrta-{self.cutoff}/{self.residue}/'
        filename = f'gibbs_{self.niter}.pkl'
        if os.path.exists(savedir):
            if os.path.exists(savedir+filename):
                os.rename(savedir+filename, savedir+filename+'.bak')
            with open(f'basicrta-{self.cutoff}/{self.residue}/gibbs_'
                      f'{self.niter}.pkl', 'w+b') as f:
                pickle.dump(self, f)
        else:
            raise OSError(f'No such directory: {savedir}')

    @staticmethod
    def load(file):
        """
        Load an instance of :class:`Gibbs`.

        :param file: Path to instance of :class:`Gibbs`
        :type file: str
        """
        from basicrta.util import get_s
        keys = ['times', 'residue', 'loc', 'ncomp', 'niter', 'g', 'burnin',
                'processed_results', 'ts', 'mcweights', 'mcrates', 't',
                's', 'cutoff', 'indicator', 'whypers', 'rhypers']
        with open(file, 'r+b') as f:
            r = pickle.load(f)

        g = Gibbs()
        for attr in keys:
            try:
                setattr(g, attr, r[f'{attr}'])
            except AttributeError:
                setattr(g, attr, None)

        if isinstance(g.residue, np.ndarray):
            g.residue = g.residue[0]

        if g.t is None:
            g.t, g.s = get_s(g.times, g.ts)

        # if len(g.processed_results) == 0:
        #     g._process_gibbs()
        return g

    def plot_tau_hist(self, scale=1, save=False):
        r"""
        Plot histogram of tau values. The figure aspect ratio is 4:3, and can be
        made larger/smaller using the `scale` argument. 

        :param scale: Increase plot size by this factor
        :type scale: float
        :param save: Save plot to file
        :type save: bool
        """
        from matplotlib.ticker import MaxNLocator
        cmap = mpl.colormaps['tab10']
        rp = self.processed_results

        imaxs = self.processed_results.indicator.max(axis=0)
        noise_inds = np.where(imaxs < self._noise_cutoff)[0]
        inds = np.delete(np.unique(rp.labels), noise_inds)
        i = rp.parameters[inds, 1].argmin()

        fig, ax = plt.subplots(1, figsize=(4*scale, 3*scale))
        ax.hist(1/rp.rates[rp.labels == i], label=f'{i}', alpha=0.5,
                color=cmap(i))
        ax.set_xlabel(r'$\tau$ [ns]')
        ax.set_ylabel('count')

        tmin = (1/rp.rates[rp.labels == i]).min()
        tmax = (1/rp.rates[rp.labels == i]).max()
        ax.set_xlim(tmin, tmax)
        ax.xaxis.set_major_locator(MaxNLocator(4))
        ax.xaxis.set_minor_locator(MaxNLocator(12))
        ax.yaxis.set_major_locator(MaxNLocator(3))
        ax.yaxis.set_minor_locator(MaxNLocator(12))
        # ax.ticklabel_format(style='sci', axis='x', scilimits=(0, 0),
        #                        useMathText=True)
        plt.tight_layout()
        if save:
            plt.savefig(f'basicrta-{self.cutoff}/{self.residue}/'
                        f'tau_hist.png',
                        bbox_inches='tight')
            plt.savefig(f'basicrta-{self.cutoff}/{self.residue}/'
                        f'tau_hist.pdf',
                        bbox_inches='tight')
        plt.show()

    def plot_hist(self, scale=1, save=False, component=None, bins=15):
        from matplotlib.ticker import MaxNLocator
        from scipy import stats
        from matplotlib.gridspec import GridSpec

        cmap = mpl.colormaps['tab10']
        rp = self.processed_results

        if component is None:
            comps = np.arange(rp.ncomp)
        elif isinstance(component, int):
            comps = [component]
        else:
            comps = component

        if self.whypers is None:
            self._prepare()

        i = comps[0]
        fig = plt.figure(figsize=(9*scale, 3*scale))
        gs = GridSpec(4, 12, figure=fig, hspace=0.2, wspace=0.2, bottom=0.28,
                      left=0.05, right=0.98, top=0.93)
        ax0 = fig.add_subplot(gs[:, :4])
        ax1 = np.array([[fig.add_subplot(gs[:-1, 4:7]),
                         fig.add_subplot(gs[:-1, 7])],
                        [fig.add_subplot(gs[-1, 4:7]),
                         fig.add_subplot(gs[-1, 7])]])
        ax2 = fig.add_subplot(gs[0, 8:]), fig.add_subplot(gs[1:, 8:])

        # plot posteriors
        [ax0.hist(rp.weights[rp.labels == i], label='posterior', alpha=0.5,
                  color=cmap(i), density=True, bins=bins) for i in comps]
        [ax1[0, 0].hist(rp.rates[rp.labels == i], label=f'{i}', alpha=0.5,
                        color=cmap(i), density=True, bins=bins) for i in comps]
        [ax1[1, 0].hist(rp.rates[rp.labels == i], label=f'{i}', alpha=0.5,
                        color=cmap(i), density=True, bins=bins) for i in comps]
        [ax2[1].hist(1/rp.rates[rp.labels == i], label=f'{i}', alpha=0.5,
                     color=cmap(i), density=True, bins=bins) for i in comps]

        # create bounds and plot priors
        wbounds = np.array([[rp.weights[rp.labels == i].min(),
                            rp.weights[rp.labels == i].max()] for i in comps])
        rbounds = np.array([[rp.rates[rp.labels == i].min(),
                            rp.rates[rp.labels == i].max()] for i in comps])
        tbounds = np.array([[(1/rp.rates[rp.labels == i]).min(),
                            (1/rp.rates[rp.labels == i]).max()] for i in comps])
        
        rx = np.linspace(0, 10, 10000)
        tx = np.linspace(0, 500, 10000)
        
        ax0.hist(rng.dirichlet(self.whypers, size=1000000).flatten(),
                 density=True, bins=20000, label='prior', alpha=0.5)
        rys = (stats.gamma(self.rhypers[0, 0], scale=1/self.rhypers[0, 1]).
               pdf(rx))
        tys = (stats.invgamma(self.rhypers[0, 0], scale=self.rhypers[0, 1]).
               pdf(tx))

        ax1[1, 0].plot(rx, rys, label=f'{i}', alpha=0.5)
        ax1[1, 0].fill_between(rx, rys, alpha=0.5)
        ax1[1, 1].plot(rx, rys, label=f'{i}', alpha=0.5)
        ax1[1, 1].fill_between(rx, rys, alpha=0.5)

        ax2[0].plot(tx, tys, label=f'{i}', alpha=0.5)
        ax2[0].fill_between(tx, tys, alpha=0.5)
        ax2[1].plot(tx, tys, label=f'{i}', alpha=0.5)
        ax2[1].fill_between(tx, tys, alpha=0.5)

        ax1[0, 0].spines['bottom'].set_visible(False)
        ax1[0, 1].spines['bottom'].set_visible(False)
        ax1[1, 0].spines['top'].set_visible(False)
        ax1[1, 1].spines['top'].set_visible(False)
        ax1[0, 0].spines['right'].set_visible(False)
        ax1[1, 0].spines['right'].set_visible(False)
        ax1[0, 1].spines['left'].set_visible(False)
        ax1[1, 1].spines['left'].set_visible(False)
        ax1[0, 0].tick_params(axis='x', labelbottom=False)
        ax1[0, 1].tick_params(axis='x', labelbottom=False)
        ax1[0, 1].tick_params(axis='y', labelleft=False)
        ax1[1, 1].tick_params(axis='y', labelleft=False)

        ax2[0].spines['bottom'].set_visible(False)
        ax2[1].spines['top'].set_visible(False)
        ax2[0].tick_params(axis='x', labelbottom=False)
        ax2[0].set_xticks([])

        d = 0.15
        kwargs = dict(marker=[(-1, -d), (1, d)], markersize=12,
                      linestyle="none", color='k', mec='k', mew=1,
                      clip_on=False)
        kwargs2 = dict(marker=[(1+d, 0), (0, 1+d)], markersize=12,
                       linestyle="none", color='k', mec='k', mew=1,
                       clip_on=False)

        ax1[0, 0].plot([0], transform=ax1[0, 0].transAxes, **kwargs)
        ax1[1, 0].plot([1], transform=ax1[1, 0].transAxes, **kwargs)
        ax1[0, 1].plot([1], [0], transform=ax1[0, 1].transAxes, **kwargs)
        ax1[1, 1].plot([1], [1], transform=ax1[1, 1].transAxes, **kwargs)
        ax1[0, 0].plot([1], [1], transform=ax1[0, 0].transAxes, **kwargs2)

        ax2[0].plot([0, 1], [0, 0], transform=ax2[0].transAxes, **kwargs)
        ax2[1].plot([0, 1], [1, 1], transform=ax2[1].transAxes, **kwargs)

        ax0.set_xlabel(r'$\pi_k$')
        ax1[1, 0].set_xlabel(r'$\lambda_k$ [ns$^{-1}$]')
        # set_shared_xlabel(ax1[1, :], label=r'$\lambda_k$ [ns$^{-1}$]')
        ax2[1].set_xlabel(r'$\tau$ [ns]')
        ax0.set_ylabel('p')
        if component is None:
            ax1[0].set_xlim(1e-4, 1)
            ax1[1].set_xlim(1e-3, 10)
            ax1[0].legend(title='component')
            ax1[1].legend(title='component')
            ax1[0].set_xscale('log')
            ax1[1].set_xscale('log')
        else:
            ax0.set_xlim(1e-5, 1e-3)

            ax1[0, 0].set_xlim(1e-4, 1e-2)
            ax1[1, 0].set_xlim(1e-4, 1e-2)
            ax1[0, 1].set_xlim(1e-2, 10)
            ax1[1, 1].set_xlim(1e-2, 10)
            ax1[0, 0].set_ylim(5, 1200)
            ax1[0, 1].set_ylim(5, 1200)
            ax1[1, 0].set_ylim(0, 5)
            ax1[1, 1].set_ylim(0, 5)

            ax2[0].set_xlim(-5, 500)
            ax2[1].set_xlim(-5, 500)
            ax2[0].set_ylim(0.02, 0.2)
            ax2[1].set_ylim(0, 0.015)

            ax0.xaxis.set_major_locator(MaxNLocator(3, min_n_ticks=3,
                                                    prune='both'))
            ax0.xaxis.set_minor_locator(MaxNLocator(12, min_n_ticks=9,
                                                    prune='both'))
            ax0.yaxis.set_major_locator(MaxNLocator(3, min_n_ticks=3,
                                                    prune='both'))
            ax0.yaxis.set_minor_locator(MaxNLocator(12, min_n_ticks=9,
                                                    prune='both'))
            ax0.ticklabel_format(style='sci', axis='both', scilimits=(0, 0),
                                 useMathText=True)

            ax1[1, 0].xaxis.set_major_locator(MaxNLocator(3, min_n_ticks=3,
                                                          prune='both'))
            ax1[1, 0].xaxis.set_minor_locator(MaxNLocator(12, min_n_ticks=9,
                                                          prune='both'))
            ax1[1, 1].xaxis.set_major_locator(MaxNLocator(3, min_n_ticks=3,
                                                          prune='both'))
            ax1[1, 1].xaxis.set_minor_locator(MaxNLocator(12, min_n_ticks=9,
                                                          prune='both'))
            ax1[0, 0].yaxis.set_major_locator(MaxNLocator(3, min_n_ticks=3,
                                                          prune='both'))
            ax1[0, 0].yaxis.set_minor_locator(MaxNLocator(12, min_n_ticks=9,
                                                          prune='both'))
            ax1[1, 0].yaxis.set_major_locator(MaxNLocator(3, min_n_ticks=3,
                                                          prune='both'))
            ax1[1, 0].yaxis.set_minor_locator(MaxNLocator(12, min_n_ticks=9,
                                                          prune='both'))
            
            ax1[0, 0].ticklabel_format(style='sci', axis='y', scilimits=(1, 1),
                                       useMathText=True)
            ax1[1, 0].ticklabel_format(style='sci', axis='y', scilimits=(1, 1),
                                       useMathText=True)
            ax1[1, 0].ticklabel_format(style='sci', axis='x', 
                                       scilimits=(-3, -3), useMathText=True)
            ax1[0, 1].ticklabel_format(style='sci', axis='x', scilimits=(0, 0),
                                       useMathText=True)
            ax1[0, 1].ticklabel_format(style='sci', axis='y', scilimits=(1, 1),
                                       useMathText=True)
            
            ax2[0].yaxis.set_major_locator(MaxNLocator(3, min_n_ticks=2,
                                                       prune='both'))
            ax2[0].yaxis.set_minor_locator(MaxNLocator(12, min_n_ticks=9,
                                                       prune='both'))
            ax2[1].yaxis.set_major_locator(MaxNLocator(3, min_n_ticks=3,
                                                       prune='both'))
            ax2[1].yaxis.set_minor_locator(MaxNLocator(15, min_n_ticks=9,
                                                       prune='both'))
            ax2[1].xaxis.set_major_locator(MaxNLocator(3, min_n_ticks=3,
                                                       prune='both'))
            ax2[1].xaxis.set_minor_locator(MaxNLocator(12, min_n_ticks=9,
                                                       prune='both'))
            
            ax2[0].ticklabel_format(style='sci', axis='x', scilimits=(0, 0),
                                    useMathText=True)
            ax2[0].ticklabel_format(style='sci', axis='y', scilimits=(-1, -1),
                                    useMathText=True)
            ax2[1].ticklabel_format(style='sci', axis='y', scilimits=(-1, -1),
                                    useMathText=True)
            ax2[1].ticklabel_format(style='sci', axis='x', scilimits=(2, 2),
                                    useMathText=True)

            ax1[0, 0].set_xticks([])
            ax1[0, 1].set_xticks([])
            ax1[0, 1].set_yticks([])
            ax1[1, 1].set_yticks([])
            handles, labels = ax0.get_legend_handles_labels()
            fig.legend(handles, labels, loc='lower center', ncols=2)
        if save:
            if component is not None:
                plt.savefig(f'basicrta-{self.cutoff}/{self.residue}/'
                            f'hist_results_{component}.png',
                            bbox_inches='tight')
                plt.savefig(f'basicrta-{self.cutoff}/{self.residue}/'
                            f'hist_results_{component}.pdf',
                            bbox_inches='tight')
            else:
                plt.savefig(f'basicrta-{self.cutoff}/{self.residue}/'
                            'hist_results.png', bbox_inches='tight')
                plt.savefig(f'basicrta-{self.cutoff}/{self.residue}/'
                            'hist_results.pdf', bbox_inches='tight')
        plt.show()

    def plot_gibbs(self, scale=1.5, sparse=1, save=False):
        cmap = mpl.colormaps['tab10']
        rp = self.processed_results

        fig, ax = plt.subplots(2, figsize=(4*scale, 3*scale), sharex=True)
        [ax[0].plot(rp.iteration[rp.labels == i][::sparse],
                    rp.weights[rp.labels == i][::sparse], '.',
                    label=f'{i}', color=cmap(i))
         for i in np.unique(rp.labels)]
        ax[0].set_yscale('log')
        ax[0].set_ylabel(r'$\pi_k$')
        [ax[1].plot(rp.iteration[rp.labels == i][::sparse],
                    rp.rates[rp.labels == i][::sparse], '.', label=f'{i}',
                    color=cmap(i)) for i in np.unique(rp.labels)]
        ax[1].set_yscale('log')
        ax[1].set_ylabel(r'\lambda_k (ns$^{-1}$)')
        ax[1].set_xlabel('sample')
        ax[0].legend(title='component')
        ax[1].legend(title='component')
        plt.tight_layout()
        if save:
            plt.savefig(f'basicrta-{self.cutoff}/{self.residue}/'
                        'plot_results.png', bbox_inches='tight')
            plt.savefig(f'basicrta-{self.cutoff}/{self.residue}/'
                        'plot_results.pdf', bbox_inches='tight')
        plt.show()

    def _estimate_params(self):
        rp = self.processed_results

        ws = [rp.weights[rp.labels == i] for i in range(rp.ncomp)]
        rs = [rp.rates[rp.labels == i] for i in range(rp.ncomp)]
        wbins = [np.exp(np.linspace(np.log(rp.weights[rp.labels == i].min()),
                                    np.log(rp.weights[rp.labels == i].max()),
                                    20))
                 for i in range(rp.ncomp)]
        rbins = [np.exp(np.linspace(np.log(rp.rates[rp.labels == i].min()),
                                    np.log(rp.rates[rp.labels == i].max()), 20))
                 for i in range(rp.ncomp)]
        wbounds = np.array([confidence_interval(d) for d in ws])
        rbounds = np.array([confidence_interval(d) for d in rs])

        whists = [np.histogram(w, bins=bins) for w, bins in zip(ws, wbins)]
        rhists = [np.histogram(r, bins=bins) for r, bins in zip(rs, rbins)]

        params = np.array([[wh[1][np.argmax(wh[0])], rh[1][np.argmax(rh[0])]]
                           for wh, rh in zip(whists, rhists)])

        setattr(rp, 'parameters', params)
        setattr(rp, 'intervals', np.array([wbounds, rbounds]))

    def estimate_tau(self):
        r"""
        Estimate the posterior maximum and confidence interval (CI) for the
        :math:`tau` distribution of the slowest process. NOTE: In the future 
        this will return an array containing :math:`tau` and CI for all
        clusters.

        :return: An array containing the posterior maximum and bounds of the
                 95% confidence interval in the format [LB, max, UB].
        :rtype: list
        """
        bintype='sqrt'
        rp = self.processed_results

        imaxs = self.processed_results.indicator.max(axis=0)
        noise_inds = np.where(imaxs < self._noise_cutoff)[0]
        inds = np.delete(np.unique(rp.labels), noise_inds)
        index = rp.parameters[inds, 1].argmin()

        taus = 1 / rp.rates[rp.labels == index]
        wts = rp.weights[rp.labels == index]
        ci = confidence_interval(taus)
        h = np.histogram(taus, bins='sqrt')
        indmax = h[0].argmax()
        val = 0.5 * (h[1][:-1][indmax] + h[1][1:][indmax])
        
        # Used for finding maximum of weight vs tau 2d distribution
        #wbins = np.histogram_bin_edges(wts, bins=bintype)
        #rbins = np.histogram_bin_edges(taus, bins=bintype)
        #vals, ws, rs = np.histogram2d(wts, taus, bins=[wbins,rbins])
        #indmax = np.unravel_index(vals.argmax(), vals.shape)
        #val = 0.5 * (rs[:-1] + rs[1:])[indmax[1]]
        return [ci[0], val, ci[1]]

    def plot_surv(self, scale=1, remove_noise=False, save=False, xlim=None,
                  ylim=(1e-6, 5), xmajor=None, xminor=None, xscale='linear',
                  yscale='log'):
        """
        Plot the survival function with the exponential mixture components where
        parameters are determined from the clustering results.

        :param scale: Modify the size of the figure by this factor
        :type scale: float
        :param remove_noise: Whether to remove noise clusters
        :type remove_noise: bool
        :param save: Whether to save the figure
        :type save: bool
        :param xlim: X-axis limits
        :type xlim: tuple
        :param ylim: Y-axis limits
        :type ylim: tuple
        :param xmajor: X-axis major tick
        :type xmajor: int
        :param xminor: X-axis minor tick
        :type xminor: int
        """
        from matplotlib.ticker import MultipleLocator, MaxNLocator

        if xmajor is None:
            maj_loc = MaxNLocator(nbins=3)
        else:
            maj_loc = MultipleLocator(xmajor)

        if xminor is None:
            min_loc = MaxNLocator(nbins=12)
        else:
            min_loc = MultipleLocator(xminor)

        cmap = mpl.colormaps['tab10']
        rp = self.processed_results
        imaxs = self.processed_results.indicator.max(axis=0)
        noise_inds = np.where(imaxs < self._noise_cutoff)[0]
        uniq_labels = np.unique(rp.labels)
        if remove_noise:
            uniq_labels = np.delete(uniq_labels, noise_inds)

        ws, rs = rp.parameters[:, 0], rp.parameters[:, 1]
        fig, ax = plt.subplots(1, figsize=(4*scale, 3*scale))
        ax.plot(self.t, self.s, '.')
        [ax.plot(self.t, ws[i]*np.exp(-rs[i]*self.t), label=f'{i}',
                 color=cmap(i)) for i in np.unique(uniq_labels)]
        ax.set_ylim(ylim)
        ax.set_xlim(xlim)
        ax.set_yscale(yscale)
        ax.set_xscale(xscale)
        ax.set_ylabel('survival function $s$')
        ax.set_xlabel(r'$t$ [ns]')
        ax.set_yticks([1, 1e-2, 1e-4])
        ax.xaxis.set_major_locator(maj_loc)
        ax.xaxis.set_minor_locator(min_loc)
        ax.legend(title='cluster')
        plt.tight_layout()
        if save:
            plt.savefig(f'basicrta-{self.cutoff}/{self.residue}/'
                        's_vs_t.png', bbox_inches='tight')
            plt.savefig(f'basicrta-{self.cutoff}/{self.residue}/'
                        's_vs_t.pdf', bbox_inches='tight')
        plt.show()


if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument('--contacts')
    parser.add_argument('--resid', type=int, default=None)
    parser.add_argument('--nproc', type=int, default=1)
    parser.add_argument('--niter', type=int, default=110000)
    parser.add_argument('--ncomp', type=int, default=15)
    args = parser.parse_args()

    contact_path = os.path.abspath(args.contacts)
    cutoff = args.contacts.split('/')[-1].strip('.pkl').split('_')[-1]

    ParallelGibbs(contact_path, nproc=args.nproc, ncomp=args.ncomp,
                  niter=args.niter).run(run_resids=args.resid)
