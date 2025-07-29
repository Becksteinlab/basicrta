"""Functions used by other modules."""

from matplotlib.ticker import MultipleLocator
from matplotlib.patches import Rectangle
from matplotlib.collections import PatchCollection
from glob import glob
from tqdm import tqdm
import ast
import multiprocessing
import os
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
import pickle
import seaborn as sns
import MDAnalysis as mda
import MDAnalysis
import warnings
#from scipy.optimize import linear_sum_assignment as lsa
#import bz2

mpl.rcParams['pdf.fonttype'] = 42


def siground(x, dec):
    return float(f'%.{dec - 1}e' % x)


def slice_trajectory(u, nslices):
    if isinstance(u, MDAnalysis.coordinates.base.FrameIteratorSliced):
        frames = np.arange(u.start, u.stop, u.step)
    elif isinstance(u, MDAnalysis.coordinates.base.FrameIteratorIndices):
        frames = u.frames
    else:
        frames = np.arange(len(u.trajectory))

    sliced_frames = np.array_split(frames, nslices)
    return sliced_frames

# CURRENTLY UNUSED
# def KL_resort(r):
#     mcweights, mcrates = r.mcweights.copy(), r.mcrates.copy()
#     indicator[:] = indicator_bak
#     Ls, niter = [L], 0
#     for j in tqdm(range(r.niter)):
#         sorts = mcweights[j].argsort()[::-1]
#         mcweights[j] = mcweights[j][sorts]
#         mcrates[j] = mcrates[j][sorts]
#
#     while niter < 10:
#         Z = np.zeros_like(z)
#         for j in tqdm(range(2000, 3000), desc='recomputing Q'):
#             tmp = mcweights[j] * mcrates[j] * np.exp(np.outer(-mcrates[j], x)).T
#             z = (tmp.T / tmp.sum(axis=1)).T
#             Z += z
#         Z = Z / 1000
#
#         for j in tqdm(range(2000, 3000), desc='resorting'):
#             tmp = mcweights[j] * mcrates[j] * np.exp(np.outer(-mcrates[j], x)).T
#             z = (tmp.T / tmp.sum(axis=1)).T
#
#             tmpsum = np.ones((ncomp, ncomp), dtype=np.float64)
#             for k in range(ncomp):
#                 tmpsum[k] = np.sum(z[:, k] * np.log(z[:, k] / Z.T), axis=1)
#
#             tmpsum[tmpsum != tmpsum] = 1e20
#             sorts = lsa(tmpsum)[1]
#             mcweights[j] = mcweights[j][sorts]
#             mcrates[j] = mcrates[j][sorts]
#         niter += 1


def tm(Prot, i):
    dif = Prot['tm{0}'.format(i)][1] - Prot['tm{0}'.format(i)][0]
    return [Prot['tm{0}'.format(i)], dif]


def confidence_interval(data, percentage=95):
    ds = np.sort(data)
    perc = np.arange(1, len(ds) + 1) / (len(ds))
    lower = (100 - percentage) / 200
    upper = (percentage + (100 - percentage) / 2) / 100

    try:
        l = ds[np.where(perc <= lower)[0][-1]]
        u = ds[np.where(perc >= upper)[0][0]]
    except IndexError:
        l = ds[0]
        u = ds[-1]

    return [l, u]


def get_bars(tau):
    maxs = tau[:, 1]
    lb, ub = tau[:, 0], tau[:, 2]
#    if lb > maxs or ub < maxs:
#        raise ValueError('bounds not valid')
    return np.array([maxs - lb, ub - maxs])


def unique_rates(ncomp, mcrates):
    mclen = len(mcrates) * 9 // 10
    means = mcrates[mclen:].mean(axis=0)
    stds = mcrates[mclen:].std(axis=0)
    lb, ub = means - stds, means + stds
    bools = np.empty([ncomp, ncomp])
    for j, mean in enumerate(means):
        for i in range(ncomp):
            bools[j, i] = ((mean < ub[i]) & (mean > lb[i]))
    sums = bools.sum(axis=0)
    deg_rts = sums[np.where(sums != 1)]
    return ncomp - len(deg_rts)


def get_s(x, ts):
    Bins = get_bins(x, ts)
    Hist = np.histogram(x, bins=Bins)
    t, s = make_surv(Hist)
    return t, s


def plot_r_vs_w(r, rrange=None, wrange=None):
    plt.close()
    plt.figure(figsize=(4, 3))
    for k in range(r.ncomp):
        plt.plot(r.mcrates[:, k], r.mcweights[:, k], label=f'{k}')
    plt.yscale('log')
    plt.xscale('log')
    if rrange:
        plt.xlim(*rrange)
    if wrange:
        plt.ylim(*wrange)
    plt.ylabel('weight')
    plt.xlabel('rate')
    plt.legend(loc='upper left')
    plt.savefig(f'{r.name}/figs/k{r.ncomp}_r_vs_w.png')
    plt.savefig(f'{r.name}/figs/k{r.ncomp}_r_vs_w.pdf')


# CURRENTLY UNUSED
# def plot_r_vs_w(weights, rates, labels, rrange=None, wrange=None):
#     plt.close()
#     plt.figure(figsize=(4, 3))
#     ncomp = len(np.unique(labels))
#     for k in range(ncomp):
#         inds = np.where(labels == k)[0]
#         plt.plot(rates[inds], weights[inds], '.', label=f'{k}')
#     plt.yscale('log')
#     plt.xscale('log')
#     if rrange:
#         plt.xlim(*rrange)
#     if wrange:
#         plt.ylim(*wrange)
#     plt.ylabel('weight')
#     plt.xlabel('rate')
#     plt.legend(loc='upper left')
#     plt.tight_layout()
#     plt.show()


def get_color(i):
    if i < 0:
        color = i
    else:
        color = i % 20
    return color


def plot_results(results, cond='ml', save=False, show=False):
    outdir = results.name
    sortinds = np.argsort([line.mean() for line in results.rates])

    weight_posts = np.array(getattr(results, 'weights'), dtype=object)[sortinds]
    rate_posts = np.array(getattr(results, 'rates'), dtype=object)[sortinds]
    w_hists = [plt.hist(post, density=True) for post in weight_posts]
    r_hists = [plt.hist(post, density=True) for post in rate_posts]
    plt.close('all')
    if cond == 'mean':
        weights = np.array([w.mean() for w in results.weights])
        weights = weights / weights.sum()
        rates = np.array([r.mean() for r in results.rates])
    elif cond == 'ml':
        mlw, mlr = [], []
        for i in range(results.ncomp):
            mlw.append(w_hists[i][1][w_hists[i][0].argmax()])
            mlr.append(r_hists[i][1][r_hists[i][0].argmax()])
        mlw = np.array(mlw)
        weights = mlw / mlw.sum()
        rates = np.array(mlr)
    else:
        raise ValueError('Only implemented for most likely (ml) and mean')

    fig, axs = plt.subplots(figsize=(4, 3))
    plt.scatter(results.t, results.s, s=15, label='data')
    plt.plot(results.t, np.inner(weights, np.exp(np.outer(results.t, -rates))),
             label='fit', color='y', ls='dashed', lw=3)
    for i in range(results.ncomp):
        plt.plot(results.t, weights[i] * np.exp(results.t * -rates[i]),
                 label=f'Comp.{i}', color=f'C{i}')
    plt.plot([], [], ' ', label=rf'$\tau$={np.round(1 / rates.min(), 1)} ns')
    plt.yscale('log')
    plt.ylim(0.8 * results.s[-2], 2)
    plt.xlim(-0.05 * results.t[-2], 1.1 * results.t[-2])
    plt.legend()
    plt.ylabel('s').set_rotation(0)
    plt.xlabel('time (ns)')
    plt.tight_layout()
    sns.despine(offset=3, ax=axs)
    if save:
        plt.savefig(f'{outdir}/figs/k{results.ncomp}-{cond}_results.png')
        plt.savefig(f'{outdir}/figs/k{results.ncomp}-{cond}_results.pdf')
    if show:
        plt.show()
    plt.close('all')


def all_post_hist(results, save=False, show=False, wlims=None, rlims=None):
    outdir = results.name
    for attr, unit in [['rates', ' (ns$^{-1}$)'], ['weights', '']]:
        Attr = getattr(results, attr)
        plt.figure(figsize=(4, 3))
        for i in range(results.ncomp):
            plt.hist(Attr[i], density=True, bins=15, label=f'comp. {i}',
                     alpha=0.5)
        plt.legend()
        plt.xlabel(f'{attr}{unit}'), plt.ylabel('p').set_rotation(0)
        plt.yscale('log'), plt.xscale('log')
        if attr == 'rates' and rlims:
            plt.xlim(rlims[0])
            plt.ylim(rlims[1])
        if attr == 'weights' and wlims:
            plt.xlim(wlims[0])
            plt.ylim(wlims[1])
        if save:
            name = f'{outdir}/figs/k{results.ncomp}-posterior_{attr}_comp-all'
            plt.savefig(f'{name}.png', bbox_inches='tight')
            plt.savefig(f'{name}.pdf', bbox_inches='tight')
        if show:
            plt.show()
        plt.close('all')


def plot_post(results, attr, comp=None, save=False, show=False):
    outdir = results.name
    Attr = getattr(results, attr)
    if attr == 'rates':
        unit = ' (ns$^{-1}$)'
    else:
        unit = ''

    if comp:
        [plt.hist(Attr[i], density=True, bins=50, label=f'comp. {i}') for i in
         comp]
        plt.legend()
        if save:
            plt.savefig(f'{outdir}/figs/k{results.ncomp}-posterior_{attr}_\
                           comps-{"-".join([str(i) for i in comp])}.png')
            plt.savefig(f'{outdir}/figs/k{results.ncomp}-posterior_{attr}_\
                           comps-{"-".join([str(i) for i in comp])}.pdf')
        if show:
            plt.show()
        plt.close('all')
    else:
        for i in range(results.ncomp):
            plt.close()
            fig, ax = plt.subplots(figsize=(4, 3))
            plt.hist(Attr[i], density=True, bins=15, label=f'comp. {i}')
            # plt.legend()
            plt.ylabel('p').set_rotation(0)
            plt.xlabel(rf'{attr[:-1]} {unit}')
            ax.xaxis.major.formatter._useMathText = True
            if save:
                plt.savefig(f'{outdir}/figs/k{results.ncomp}-posterior_{attr}_'
                            'comp-{i}.png', bbox_inches='tight')
                plt.savefig(f'{outdir}/figs/k{results.ncomp}-posterior_{attr}_'
                            'comp-{i}.pdf', bbox_inches='tight')
            if show:
                plt.show()


def plot_trace(results, attr, comp=None, xrange=None, yrange=None, save=False,
               show=False):
    outdir = results.name
    if attr == 'weights':
        tmp = getattr(results, 'mcweights')
    elif attr == 'rates':
        tmp = getattr(results, 'mcrates')
    if not comp:
        plt.figure(figsize=(4, 3))
        for j in range(results.ncomp):
            plt.plot(range(tmp.shape[0]), tmp[:, j], label=f'Comp. {j}')
        plt.xlabel('iteration')
        plt.ylabel(f'{attr}')
        plt.legend()
        if xrange is not None:
            plt.xlim(xrange[0], xrange[1])
        if yrange is not None:
            plt.ylim(yrange[0], yrange[1])
        if save:
            plt.savefig(f'{outdir}/figs/k{results.ncomp}-trace_{attr}.png')
            plt.savefig(f'{outdir}/figs/k{results.ncomp}-trace_{attr}.pdf')
    if comp:
        plt.figure(figsize=(4, 3))
        for i in comp:
            plt.plot(range(tmp.shape[0]), tmp[:, i], label=f'Comp. {i}')
            plt.xlabel('iteration')
            plt.ylabel(f'{attr}')
            plt.legend()
        if xrange is not None:
            plt.xlim(xrange[0], xrange[1])
        if yrange is not None:
            plt.ylim(yrange[0], yrange[1])
        if save:
            plt.savefig(f'{outdir}/figs/k{results.ncomp}-trace_{attr}_comps-\
                          {"-".join([str(i) for i in comp])}.png')
            plt.savefig(f'{outdir}/figs/k{results.ncomp}-trace_{attr}_comps-\
                          {"-".join([str(i) for i in comp])}.pdf')
    if show:
        plt.show()
    plt.close('all')

# CURRENTLY UNUSED
# def collect_results(ncomp=None):
#     """returns (residues, tslow, stds)
#     """
#     dirs = np.array(glob('?[0-9]*'))
#     sorted_inds = np.array([int(adir[1:]) for adir in dirs]).argsort()
#     dirs = dirs[sorted_inds]
#     t_slow = np.zeros(len(dirs))
#     sd = np.zeros((len(dirs), 2))
#     residues = np.empty((len(dirs)), dtype=object)
#     indicators = []
#     for i, adir in enumerate(tqdm(dirs, desc='Collecting results')):
#         residues[i] = adir
#         try:
#             tmp_res = pickle.load(
#                 bz2.BZ2File(f'{adir}/results_20000.pkl.bz2', 'rb'))
#             tmp_res, rpinds = process_gibbs(tmp_res)
#         #    with open(f'{adir}/processed_results_10000.pkl', 'rb') as f:
#         #        tmp_res = pickle.load(f)
#         #    results = glob(f'{adir}/*results.pkl')
#         #    results.sort()
#         #    if ncomp and ncomp-1<=len(results):
#         #        max_comp_res = results[ncomp-2]
#         #    else:
#         #        max_comp_res = results[-1]
#         except FileNotFoundError:
#             t_slow[i] = 0
#             continue
#         # with open(max_comp_res, 'rb') as W:
#         #    tmp_res = pickle.load(W)
#
#         means = np.array([(1 / post).mean() for post in tmp_res.rates.T])
#         if len(means) == 0:
#             continue
#         ind = np.where(means == means.max())[0][0]
#         t_slow[i] = means[ind]
#         sd[i] = get_bars(1 / tmp_res.rates.T[ind])
#         indicators.append(
#             (tmp_res.indicator.T / tmp_res.indicator.sum(axis=1)).T)
#     return residues, t_slow, sd.T, indicators


def collect_n_plot(resids, comps):
    dirs = np.array(glob('?[0-9]*'))
    tmpresids = np.array([int(adir[1:]) for adir in dirs])
    sorted_inds = tmpresids.argsort()
    tmpresids.sort()
    dirs = dirs[sorted_inds]
    idinds = np.array([np.where(tmpresids == resid)[0][0] for resid in resids])
    dirs = dirs[idinds]

    for i, adir in enumerate(tqdm(dirs, desc='Collecting results')):
        results = glob(f'{adir}/*results.pkl')
        results.sort()
        # max_comp_res = results[-1]
        for res in results:
            with open(res, 'rb') as W:
                tmp_res = pickle.load(W)

            make_residue_plots(tmp_res, comps)
            all_post_hist(tmp_res, save=True, rlims=[[1e-3, 10], [1e-2, 1e3]],
                          wlims=[[1e-4, 1.1], [1e-1, 1e4]])
            plot_r_vs_w(tmp_res, rrange=[1e-3, 10], wrange=[1e-4, 5])


def make_residue_plots(results, comps=None, show=False):
    r = results

    if not os.path.exists(f'{r.name}/figs'):
        os.mkdir(f'{r.name}/figs/')

    plot_results(r, cond='mean', save=True, show=show)
    plot_results(r, cond='ml', save=True, show=show)
    plot_post(r, 'weights', comp=comps, save=True, show=show)
    plot_post(r, 'rates', comp=comps, save=True, show=show)
    plot_trace(r, 'weights', comp=comps, save=True, show=show,
               yrange=[-0.1, 1.1])
    plot_trace(r, 'rates', comp=comps, save=True, show=show, yrange=[-0.1, 6])


def plot_protein(residues, t_slow, bars, prot=None, label_cutoff=3, ylim=None,
                 major_tick=None, minor_tick=None, scale=1):
    try:
        with open('tm_dict.txt', 'r') as f:
            contents = f.read()
            prots = ast.literal_eval(contents)
    except FileNotFoundError:
        warnings.warn("tm_dict.txt not found, TM bars will not be drawn in " 
                      r"$\tau$ vs resid plot")
        prot = None

    if not os.path.exists('figs'):
        os.mkdir('figs')

    height, width = 3 * scale, 4 * scale
    fig, axs = plt.subplots(2, 1, figsize=(width, height), sharex=True)
    axs[0].tick_params(axis='both', which='major', labelsize=10)
    axs[1].tick_params(axis='both', which='major', labelsize=10)
    if prot is not None:
        p = [Rectangle((tm(prots[prot]['helices'], i + 1)[0][0], 0),
                       tm(prots[prot]['helices'], i + 1)[1], 1, fill=True) 
             for i in range(7)]
        patches = PatchCollection(p)
        patches.set_color('C0')
        axs[1].add_collection(patches)
    
    resids = np.array([int(res[1:]) for res in residues])
    max_inds = np.where(t_slow > label_cutoff * t_slow.mean())
    axs[0].plot(resids, t_slow, '.', color='C0')
    axs[0].errorbar(resids, t_slow, yerr=bars, fmt='none', color='C0',
                    alpha=0.5)
    [axs[0].text(resids[ind], t_slow[ind], residues[ind]) for ind in
     max_inds[0]]
    axs[0].set_ylabel(r'$\tau$ [ns]')
    axs[1].set_xlabel(r'residue')
    axs[0].get_xaxis().set_visible(False)
    axs[1].get_yaxis().set_visible(False)
    axs[1].xaxis.set_major_locator(MultipleLocator(100))
    if major_tick is not None:
        axs[0].yaxis.set_major_locator(MultipleLocator(major_tick))
    if minor_tick is not None:
        axs[0].yaxis.set_minor_locator(MultipleLocator(minor_tick))
    if ylim is not None:
        axs[0].set_ylim(ylim)
    axs[1].set_aspect(7)
    axs[0].margins(x=0)

    plt.subplots_adjust(hspace=-0.45, top=0.92)
    sns.despine(offset=10, ax=axs[0], bottom=True)
    sns.despine(ax=axs[1], top=True, bottom=False, left=True)
    plt.savefig('figs/t_slow.png', bbox_inches='tight')
    plt.savefig('figs/t_slow.pdf', bbox_inches='tight')


# def plot_frame_comp(indicators, trajtimes):
#     if not os.path.exists('figs'):
#         os.mkdir('figs')
#
#     plt.scatter(np.concatenate([*trajtimes]), indicators, s=2)
#     plt.ylabel('Component')
#     plt.xlabel('Frame')
#     sns.despine(offset=10)
#     plt.tight_layout()
#     plt.savefig('figs/frame_comp.png')
#     plt.savefig('figs/frame_comp.pdf')
#  ##  plt.show()


# def run(gib):
#     gib.run()


def run_residue(residue, time, proc, ncomp, niter, cutoff, g, from_combined=False):
    from basicrta.gibbs import Gibbs
    x = np.array(time)
    if len(x) != 0:
        try:
            proc = int(multiprocessing.current_process().name.split('-')[-1])
        except ValueError:
            proc = 1

    gib = Gibbs(x, residue, proc, ncomp=ncomp, niter=niter, cutoff=cutoff)
    gib._from_combined_contacts = from_combined
    gib.run()


def check_results(residues, times, ts):
    if not os.path.exists('result_check'):
        os.mkdir('result_check')
    for time, residue in zip(times, residues):
        if os.path.exists(residue):
            kmax = glob(f'{residue}/K*_results.pkl')[-1].split('/')[-1]. \
                       split('/')[-1].split('_')[0][1:]
            os.popen(f'cp {residue}/figs/k{kmax}-mean_results.png result_check/\
                       {residue}-k{kmax}-results.png')
        else:
            t, s = get_s(np.array(time), ts)
            plt.scatter(t, s, label='data')
            plt.ylabel('s')
            plt.xlabel('t (ns)')
            plt.legend()
            plt.title('Results unavailable')
            plt.savefig(f'result_check/{residue}-s-vs-t.png')
            plt.close('all')


def get_dec(ts):
    if len(str(float(ts)).split('.')[1].rstrip('0')) == 0:
        dec = -len(str(ts)) + 1
    else:
        dec = len(str(float(ts)).split('.')[1].rstrip('0'))
    return dec


def get_start_stop_frames(simtime, timelen, ts):
    dec = get_dec(ts)
    framec = (np.round(timelen, dec) / ts).astype(int)
    frame = (np.round(simtime, dec) / ts).astype(int)
    return frame, frame + framec - 1


def get_write_frames(u, time, trajtime, lipind, comp):
    dt, comp = u.trajectory.ts.dt / 1000, comp - 2  # nanoseconds
    bframes, eframes = get_start_stop_frames(trajtime, time, dt)
    sortinds = bframes.argsort()
    bframes.sort()
    eframes, lind = eframes[sortinds], lipind[sortinds]
    tmp = [np.arange(b, e) for b, e in zip(bframes, eframes)]
    tmpL = [np.ones_like(np.arange(b, e)) * l for b, e, l in
            zip(bframes, eframes, lind)]
    write_frames, write_Linds = np.concatenate([*tmp]), np.concatenate(
        [*tmpL]).astype(int)
    return write_frames, write_Linds


def write_trajs(u, time, trajtime, indicator, residue, lipind, step):
    try:
        proc = int(multiprocessing.current_process().name[-1])
    except ValueError:
        proc = 1

    prot, chol = u.select_atoms('protein'), u.select_atoms('resname CHOL')
    # dt = u.trajectory.ts.dt/1000 #nanoseconds
    inds = np.array(
        [np.where(indicator.argmax(axis=0) == i)[0] for i in range(8)],
        dtype=object)
    lens = np.array([len(ind) for ind in inds])
    for comp in np.where(lens != 0)[0]:
        write_frames, write_Linds = get_write_frames(u, time, trajtime, lipind,
                                                     comp + 2)
        if len(write_frames) > step:
            write_frames = write_frames[::step]
            write_Linds = write_Linds[::step]
        with mda.Writer(f"{residue}/comp{comp}_traj.xtc",
                        len((prot + chol.residues[0].atoms).atoms)) as W:
            for i, ts in tqdm(enumerate(u.trajectory[write_frames]),
                              desc=f"{residue}-comp{comp}", position=proc,
                              leave=False, total=len(write_frames)):
                ag = prot + chol.residues[write_Linds[i]].atoms
                W.write(ag)


def plot_hists(timelens, indicators, residues):
    for timelen, indicator, residue in tqdm(zip(timelens, indicators, residues),
                                            total=len(timelens),
                                            desc='ploting hists'):
        ncomps = indicator[:, 0].shape[0]

        plt.close()
        for i in range(ncomps):
            h, edges = np.histogram(timelen, density=True, bins=50,
                                    weights=indicator[i])
            m = 0.5 * (edges[1:] + edges[:-1])
            plt.plot(m, h, '.', label=i, alpha=0.5)
        plt.ylabel('p')
        plt.xlabel('time (ns)')
        plt.xscale('log')
        plt.yscale('log')
        plt.ylim(1e-6, 1)
        sns.despine(offset=5)
        plt.legend()
        plt.savefig(f'result_check/{residue}_hists_{ncomps}.png')
        plt.savefig(f'result_check/{residue}_hists_{ncomps}.pdf')


def get_remaining_residue_inds(residues, invert=True):
    dirs = np.array(glob('?[0-9]*'))
    resids = np.array([int(res[1:]) for res in residues])
    rem_residues = np.setdiff1d(residues, dirs)
    rem_resids = np.array([int(res[1:]) for res in rem_residues])
    rem_inds = np.in1d(resids, rem_resids, invert=invert)
    return rem_inds


def simulate_hn(n, weights, rates):
    n = int(n)
    x = np.zeros(n)
    p = np.random.rand(n)

    tmpw = np.concatenate(([0], np.cumsum(weights)))
    for i in range(len(weights)):
        x[(p > tmpw[i]) & (p <= tmpw[i + 1])] = \
            -np.log(
                np.random.rand(len(p[(p > tmpw[i]) & (p <= tmpw[i + 1])]))) / \
            rates[i]
    x.sort()
    return x


def make_surv(ahist):
    y = ahist[0][ahist[0] != 0]
    tmpbin = ahist[1][:-1]
    t = tmpbin[ahist[0] != 0]
    t = np.insert(t, 0, 0)
    y = np.cumsum(y)
    y = np.insert(y, 0, 0)
    y = y / y[-1]
    s = 1 - y
    return t, s


def expand_times(contacts):
    a = contacts
    prots = np.unique(a[:, 0])
    lips = np.unique(a[:, 1])

    restimes = []
    Ns = []
    for i in tqdm(prots, desc='expanding times'):
        liptimes = []
        lipNs = []
        for j in lips:
            tmp = a[(a[:, 0] == i) & (a[:, 1] == j)]
            liptimes.append(np.round(tmp[:, 2], 1))
            lipNs.append(tmp[:, 3])
        restimes.append(liptimes)
        Ns.append(lipNs)
    times = np.asarray(restimes)
    Ns = np.asarray(Ns)

    alltimes = []
    for res in tqdm(range(times.shape[0])):
        restimes = []
        for lip in range(times.shape[1]):
            for i in range(times[res, lip].shape[0]):
                [restimes.append(j) for j in [times[res, lip][i]] *
                 Ns[res, lip][i].astype(int)]
        alltimes.append(restimes)
    return np.asarray(alltimes)


def get_bins(x, ts):
    if isinstance(x, list):
        x = np.asarray(x)
    elif isinstance(x, np.ndarray):
        pass
    else:
        raise TypeError('Input should be a list or array')
    return np.arange(1, int(x.max() // ts) + 3) * ts


def extract_data(gibbs):
    from scipy import stats

    burnin_ind = gibbs.burnin // gibbs.g
    data_len = len(gibbs.times)
    wcutoff = 10 / data_len

    weights, rates = gibbs.mcweights[burnin_ind:], gibbs.mcrates[burnin_ind:]
    lens = np.array([len(row[row > wcutoff]) for row in weights])
    lmin, lmode, lmax = lens.min(), stats.mode(lens).mode, lens.max()
    train_param = lmode

    train_inds = np.where(lens == train_param)[0]
    train_weights = (weights[train_inds][weights[train_inds] > wcutoff].
                     reshape(-1, train_param))
    train_rates = (rates[train_inds][weights[train_inds] > wcutoff].
                   reshape(-1, train_param))

    inds = np.where(weights > wcutoff)
    aweights, arates = weights[inds], rates[inds]
    rcutoff = arates.min()
    data = np.stack((aweights, arates), axis=1)

    # tweights, trates = train_weights.flatten(), train_rates.flatten()
    # train_data = np.stack((tweights, trates), axis=1)

    # tmpw, tmpr = np.delete(weights, train_inds), np.delete(rates, train_inds)
    # pweights, prates = tmpw[tmpw > wcutoff], tmpr[tmpw > wcutoff]
    # predict_data = np.stack((pweights, prates), axis=1)
    return data, train_inds

def mixture_and_plot(gibbs, scale=2, sparse=1, remove_noise=False, wlim=None,
                     rlim=None, show=True, **kwargs):
    from scipy import stats
    from matplotlib.ticker import MaxNLocator

    burnin_ind = gibbs.burnin // gibbs.g
    data_len = len(gibbs.times)
    wcutoff = 10 / data_len
    if wlim is not None:
        wmin, wmax = wlim[0], wlim[1]
    else:
        wmin, wmax = wcutoff, 2

    weights = gibbs.mcweights[burnin_ind::gibbs.gskip] 
    rates = gibbs.mcrates[burnin_ind::gibbs.gskip]
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

    rcutoff = arates.min()
    if rlim is not None:
        rmin, rmax = rlim[0], rlim[1]
    else:
        rmin, rmax = rcutoff, 10

    tweights, trates = train_weights.flatten(), train_rates.flatten()
    train_data = np.stack((tweights, trates), axis=1)

    tmpw, tmpr = np.delete(weights, train_inds), np.delete(rates, train_inds)
    pweights, prates = tmpw[tmpw > wcutoff], tmpr[tmpw > wcutoff]
    predict_data = np.stack((pweights, prates), axis=1)

    train_inds = np.array([np.where(data == val)[0] for val in train_data])
    predict_inds = np.array([np.where(data == val)[0] for val in predict_data])

    all_labels = gibbs.processed_results.labels
    uniq_labels = np.unique(all_labels)
    leg_labels = np.array([' ' for _ in uniq_labels])
    labels = all_labels[train_inds]
    predict_labels = all_labels[predict_inds]

    imaxs = gibbs.processed_results.indicator.max(axis=0)
    noise_inds = np.where(imaxs < gibbs._noise_cutoff)[0]

    means = np.array([arates[all_labels == i].mean() for i in uniq_labels])
    vsorts = means[np.delete(uniq_labels, noise_inds)].argsort()[::-1]
    nsorts = means[noise_inds].argsort()[::-1]
    presorts = np.concatenate([np.delete(uniq_labels, noise_inds)[vsorts],
                               noise_inds[nsorts]])
    sorts = np.array([np.where(presorts == i)[0][0] for i in uniq_labels])

    labels = sorts[labels]
    predict_labels = sorts[predict_labels]
    all_labels = sorts[all_labels]
    uniq_vlabels = uniq_labels[:len(vsorts)]
    if remove_noise:
        uniq_labels = uniq_vlabels

    tinds = [np.where(labels == i)[0] for i in uniq_labels]
    pinds = [np.where(predict_labels == i)[0] for i in uniq_labels]

    train_data_inds = np.array([np.where(data == col)[0][0] for col in
                                train_data])
    predict_data_inds = np.array([np.where(data == col)[0][0] for col in
                                  predict_data])

    cmap = mpl.colormaps['tab10']
    cmap.set_under()

    figa, axa = plt.subplots(2, 2, figsize=(4 * scale, 3 * scale))
    figt, axt = plt.subplots(2, 2, figsize=(4 * scale, 3 * scale))
    figp, axp = plt.subplots(2, 2, figsize=(4 * scale, 3 * scale))

    # create histogram
    fig1a, ax1a = plt.subplots(1, figsize=(4, 3))
    fig1t, ax1t = plt.subplots(1, figsize=(4, 3))
    fig1p, ax1p = plt.subplots(1, figsize=(4, 3))
    for i in uniq_labels[::-1]:
        # bins = np.exp(np.linspace(np.log(trates[tinds[i]].min()),
        #                           np.log(trates[tinds[i]].max()), 50))
        bins = np.linspace(trates[tinds[i]].min(), trates[tinds[i]].max(), 50)
        axa[0, 0].hist(prates[pinds[i]], bins=bins, label=leg_labels[i],
                       color=cmap(get_color(i)), zorder=1, alpha=0.5)
        axa[0, 0].hist(trates[tinds[i]], bins=bins, label=leg_labels[i],
                       color=cmap(get_color(i)), zorder=2, alpha=0.5,
                       edgecolor='k')
        axt[0, 0].hist(trates[tinds[i]], bins=bins, label=leg_labels[i],
                       color=cmap(get_color(i)), zorder=2, alpha=0.5,
                       edgecolor='k')
        axp[0, 0].hist(prates[pinds[i]], bins=bins, label=leg_labels[i],
                       color=cmap(get_color(i)), zorder=1, alpha=0.5)

        ax1a.hist(prates[pinds[i]], bins=bins, label=leg_labels[i],
                  color=cmap(get_color(i)), zorder=1, alpha=0.5)
        ax1a.hist(trates[tinds[i]], bins=bins, label=leg_labels[i],
                  color=cmap(get_color(i)), zorder=2, alpha=0.5,
                  edgecolor='k')
        ax1p.hist(prates[pinds[i]], bins=bins, label=leg_labels[i],
                  color=cmap(get_color(i)), zorder=1, alpha=0.5)
        ax1t.hist(trates[tinds[i]], bins=bins, label=leg_labels[i],
                  color=cmap(get_color(i)), zorder=2, alpha=0.5,
                  edgecolor='k')

    # for combined plot
    locatora = MaxNLocator(prune='both', nbins=5, min_n_ticks=3)
    locatort = MaxNLocator(prune='both', nbins=5, min_n_ticks=3)
    locatorp = MaxNLocator(prune='both', nbins=5, min_n_ticks=3)

    axa[0, 0].set_xscale('log')
    axa[0, 0].set_xlabel(r'$\lambda_k$ [ns$^{-1}$]')
    axa[0, 0].set_ylabel('count')
    axa[0, 0].set_xlim(rmin, rmax)
    axa[0, 0].set_ylim(bottom=wmin)
    axa[0, 0].yaxis.set_major_locator(locatora)

    axt[0, 0].set_xscale('log')
    axt[0, 0].set_xlabel(r'$\lambda_k$ [ns$^{-1}$]')
    axt[0, 0].set_ylabel('count')
    axt[0, 0].set_xlim(rmin, rmax)
    axt[0, 0].set_ylim(bottom=wmin)
    axt[0, 0].yaxis.set_major_locator(locatort)

    axp[0, 0].set_xscale('log')
    axp[0, 0].set_xlabel(r'$\lambda_k$ [ns$^{-1}$]')
    axp[0, 0].set_ylabel('count')
    axp[0, 0].set_xlim(rmin, rmax)
    axp[0, 0].set_ylim(bottom=wmin)
    axp[0, 0].yaxis.set_major_locator(locatorp)

    # for individual plot
    ax1a.set_xscale('log')
    ax1a.set_xlabel(r'$\lambda_k$ [ns$^{-1}$]')
    ax1a.set_ylabel('count')
    ax1a.set_xlim(rmin, rmax)
    ax1a.set_ylim(bottom=wmin)
    ax1a.yaxis.set_major_locator(locatora)

    ax1t.set_xscale('log')
    ax1t.set_xlabel(r'$\lambda_k$ [ns$^{-1}$]')
    ax1t.set_ylabel('count')
    ax1t.set_xlim(rmin, rmax)
    ax1t.set_ylim(bottom=wmin)
    ax1t.yaxis.set_major_locator(locatort)

    ax1p.set_xscale('log')
    ax1p.set_xlabel(r'$\lambda_k$ [ns$^{-1}$]')
    ax1p.set_ylabel('count')
    ax1p.set_xlim(rmin, rmax)
    ax1p.set_ylim(bottom=wmin)
    ax1p.yaxis.set_major_locator(locatorp)

    for suffix in ['png', 'pdf']:
        basename = f"basicrta-{gibbs.cutoff}/{gibbs.residue}/result_hist"
        if remove_noise:
            fig1a.savefig(f"{basename}_all_noiserm.{suffix}",
                          bbox_inches='tight')
            fig1t.savefig(f"{basename}_train_noiserm.{suffix}",
                          bbox_inches='tight')
            fig1p.savefig(f"{basename}_validate_noiserm.{suffix}",
                          bbox_inches='tight')
        else:
            fig1a.savefig(f"{basename}_all.{suffix}",
                          bbox_inches='tight')
            fig1t.savefig(f"{basename}_train.{suffix}",
                          bbox_inches='tight')
            fig1p.savefig(f"{basename}_validate.{suffix}",
                          bbox_inches='tight')
    for fig in [fig1a, fig1t, fig1p]:
        plt.close(fig=fig)

    # create weight, rate vs sample plots
    fig1a, ax1a = plt.subplots(1, figsize=(4, 3))
    fig1t, ax1t = plt.subplots(1, figsize=(4, 3))
    fig1p, ax1p = plt.subplots(1, figsize=(4, 3))
    fig2a, ax2a = plt.subplots(1, figsize=(4, 3))
    fig2t, ax2t = plt.subplots(1, figsize=(4, 3))
    fig2p, ax2p = plt.subplots(1, figsize=(4, 3))

    row, col = gibbs.mcweights[burnin_ind:].shape
    iter_arr = np.mgrid[:row, :col][0]
    iters = iter_arr[inds]
    titer, piter = iters[train_data_inds], iters[predict_data_inds]

    for i in uniq_labels[::-1]:
        axa[0, 1].plot(piter[pinds[i]], pweights[pinds[i]][::sparse], '.',
                       label=leg_labels[i], color=cmap(get_color(i)), zorder=1)
        axa[1, 1].plot(piter[pinds[i]], prates[pinds[i]][::sparse], '.',
                       label=leg_labels[i], color=cmap(get_color(i)), zorder=1)
        axa[0, 1].plot(titer[tinds[i]], tweights[tinds[i]][::sparse], '.',
                       label=leg_labels[i], color=cmap(get_color(i)), zorder=2,
                       alpha=0.5, mec='k', mew=0.5)
        axa[1, 1].plot(titer[tinds[i]], trates[tinds[i]][::sparse], '.',
                       label=leg_labels[i], color=cmap(get_color(i)), zorder=2,
                       alpha=0.5, mec='k', mew=0.5)
        axt[0, 1].plot(titer[tinds[i]], tweights[tinds[i]][::sparse], '.',
                       label=leg_labels[i], color=cmap(get_color(i)), zorder=2,
                       alpha=0.9, mec='k', mew=0.5)
        axt[1, 1].plot(titer[tinds[i]], trates[tinds[i]][::sparse], '.',
                       label=leg_labels[i], color=cmap(get_color(i)), zorder=2,
                       alpha=0.9, mec='k', mew=0.5)
        axp[0, 1].plot(piter[pinds[i]], pweights[pinds[i]][::sparse], '.',
                       alpha=0.5, label=leg_labels[i], color=cmap(get_color(i)),
                       zorder=1)
        axp[1, 1].plot(piter[pinds[i]], prates[pinds[i]][::sparse], '.',
                       alpha=0.5, label=leg_labels[i], color=cmap(get_color(i)),
                       zorder=1)

        ax1a.plot(piter[pinds[i]], pweights[pinds[i]][::sparse], '.',
                  label=leg_labels[i], color=cmap(get_color(i)), zorder=1)
        ax1a.plot(titer[tinds[i]], tweights[tinds[i]][::sparse], '.',
                  label=leg_labels[i], color=cmap(get_color(i)), zorder=2,
                  alpha=0.5, mec='k', mew=0.5)
        ax2a.plot(piter[pinds[i]], prates[pinds[i]][::sparse], '.',
                  label=leg_labels[i], color=cmap(get_color(i)), zorder=1)
        ax2a.plot(titer[tinds[i]], trates[tinds[i]][::sparse], '.',
                  label=leg_labels[i], color=cmap(get_color(i)), zorder=2,
                  alpha=0.5, mec='k', mew=0.5)

        ax1t.plot(titer[tinds[i]], tweights[tinds[i]][::sparse], '.',
                  label=leg_labels[i], color=cmap(get_color(i)), zorder=2,
                  alpha=0.9, mec='k', mew=0.5)
        ax1p.plot(piter[pinds[i]], pweights[pinds[i]][::sparse], '.',
                  alpha=0.7, label=leg_labels[i], color=cmap(get_color(i)),
                  zorder=1)

        ax2t.plot(titer[tinds[i]], trates[tinds[i]][::sparse], '.',
                  label=leg_labels[i], color=cmap(get_color(i)), zorder=2,
                  alpha=0.9, mec='k', mew=0.5)
        ax2p.plot(piter[pinds[i]], prates[pinds[i]][::sparse], '.',
                  alpha=0.7, label=leg_labels[i], color=cmap(get_color(i)),
                  zorder=1)

    # for combined plots
    axa[0, 1].set_yscale('log')
    axa[0, 1].set_ylabel(r'$\pi_k$')
    axa[1, 1].set_yscale('log')
    axa[1, 1].set_ylabel(r'$\lambda_k$ [ns$^{-1}$]')
    axa[1, 1].set_xlabel('sample')
    axa[0, 1].set_xlabel('sample')
    axa[0, 1].set_ylim(wmin, wmax)
    axa[1, 1].set_xlabel('sample')
    axa[1, 1].set_ylim(rmin, rmax)

    axt[0, 1].set_yscale('log')
    axt[0, 1].set_ylabel(r'$\pi_k$')
    axt[1, 1].set_yscale('log')
    axt[1, 1].set_ylabel(r'$\lambda_k$ [ns$^{-1}$]')
    axt[1, 1].set_xlabel('sample')
    axt[0, 1].set_xlabel('sample')
    axt[0, 1].set_ylim(wmin, wmax)
    axt[1, 1].set_xlabel('sample')
    axt[1, 1].set_ylim(rmin, rmax)

    axp[0, 1].set_yscale('log')
    axp[0, 1].set_ylabel(r'$\pi_k$')
    axp[1, 1].set_yscale('log')
    axp[1, 1].set_ylabel(r'$\lambda_k$ [ns$^{-1}$]')
    axp[1, 1].set_xlabel('sample')
    axp[0, 1].set_xlabel('sample')
    axp[0, 1].set_ylim(wmin, wmax)
    axp[1, 1].set_xlabel('sample')
    axp[1, 1].set_ylim(rmin, rmax)

    # for individual plots
    ax1a.set_yscale('log')
    ax1a.set_ylabel(r'$\pi_k$')
    ax1a.set_ylim(wmin, wmax)
    ax1a.set_xlabel('sample')

    ax2a.set_yscale('log')
    ax2a.set_ylabel(r'$\lambda_k$ [ns$^{-1}$]')
    ax2a.set_ylim(rmin, rmax)
    ax2a.set_xlabel('sample')

    ax1t.set_yscale('log')
    ax1t.set_ylabel(r'$\pi_k$')
    ax1t.set_ylim(wmin, wmax)
    ax1t.set_xlabel('sample')

    ax2t.set_yscale('log')
    ax2t.set_ylabel(r'$\lambda_k$ [ns$^{-1}$]')
    ax2t.set_ylim(rmin, rmax)
    ax2t.set_xlabel('sample')

    ax1p.set_yscale('log')
    ax1p.set_ylabel(r'$\pi_k$')
    ax1p.set_ylim(wmin, wmax)
    ax1p.set_xlabel('sample')

    ax2p.set_yscale('log')
    ax2p.set_ylabel(r'$\lambda_k$ [ns$^{-1}$]')
    ax2p.set_ylim(rmin, rmax)
    ax2p.set_xlabel('sample')

    for suffix in ['png', 'pdf']:
        wbasename = f"basicrta-{gibbs.cutoff}/{gibbs.residue}/weight_results"
        rbasename = f"basicrta-{gibbs.cutoff}/{gibbs.residue}/rate_results"
        if remove_noise:
            fig1a.savefig(f"{wbasename}_all_noiserm.{suffix}",
                          bbox_inches='tight')
            fig1t.savefig(f"{wbasename}_train_noiserm.{suffix}",
                          bbox_inches='tight')
            fig1p.savefig(f"{wbasename}_validate_noiserm.{suffix}",
                          bbox_inches='tight')
            fig2a.savefig(f"{rbasename}_all_noiserm.{suffix}",
                          bbox_inches='tight')
            fig2t.savefig(f"{rbasename}_train_noiserm.{suffix}",
                          bbox_inches='tight')
            fig2p.savefig(f"{rbasename}_validate_noiserm.{suffix}",
                          bbox_inches='tight')
        else:
            fig1a.savefig(f"{wbasename}_all.{suffix}",
                          bbox_inches='tight')
            fig1t.savefig(f"{wbasename}_train.{suffix}",
                          bbox_inches='tight')
            fig1p.savefig(f"{wbasename}_validate.{suffix}",
                          bbox_inches='tight')
            fig2a.savefig(f"{rbasename}_all.{suffix}",
                          bbox_inches='tight')
            fig2t.savefig(f"{rbasename}_train.{suffix}",
                          bbox_inches='tight')
            fig2p.savefig(f"{rbasename}_validate.{suffix}",
                          bbox_inches='tight')

    for fig in [fig1a, fig1t, fig1p, fig2a, fig2t, fig2p]:
        plt.close(fig=fig)

    # create weight vs rate plot
    fig1a, ax1a = plt.subplots(1, figsize=(4, 3))
    fig1t, ax1t = plt.subplots(1, figsize=(4, 3))
    fig1p, ax1p = plt.subplots(1, figsize=(4, 3))
    for i in uniq_labels[::-1]:
        axa[1, 0].plot(prates[pinds[i]], pweights[pinds[i]], '.',
                       alpha=0.7, label=leg_labels[i],
                       color=cmap(get_color(i)), zorder=1)
        axa[1, 0].plot(trates[tinds[i]], tweights[tinds[i]], '.',
                       label=leg_labels[i],
                       color=cmap(get_color(i)), zorder=2, alpha=0.5,
                       mec='k', mew=0.5)
        axt[1, 0].plot(trates[tinds[i]], tweights[tinds[i]], '.',
                       label=leg_labels[i],
                       color=cmap(get_color(i)), zorder=2, alpha=0.5,
                       mec='k', mew=0.5)
        axp[1, 0].plot(prates[pinds[i]], pweights[pinds[i]], '.',
                       alpha=0.7, label=leg_labels[i],
                       color=cmap(get_color(i)), zorder=1)

        ax1a.plot(prates[pinds[i]], pweights[pinds[i]], '.', alpha=0.7,
                  label=leg_labels[i], color=cmap(get_color(i)), zorder=1)
        ax1a.plot(trates[tinds[i]], tweights[tinds[i]], '.',
                  label=leg_labels[i], color=cmap(get_color(i)), zorder=2,
                  alpha=0.5, mec='k', mew=0.5)
        ax1t.plot(trates[tinds[i]], tweights[tinds[i]], '.',
                  label=leg_labels[i], color=cmap(get_color(i)), zorder=2,
                  alpha=0.5, mec='k', mew=0.5)
        ax1p.plot(prates[pinds[i]], pweights[pinds[i]], '.', alpha=0.7,
                  label=leg_labels[i], color=cmap(get_color(i)), zorder=1)

    # for combined plots
    axa[1, 0].set_yscale('log')
    axa[1, 0].set_xscale('log')
    axa[1, 0].set_ylabel(r'$\pi_k$')
    axa[1, 0].set_xlabel(r'$\lambda_k$ [ns$^{-1}$]')
    axa[1, 0].set_xlim(rmin, rmax)
    axa[1, 0].set_ylim(wmin, wmax)

    axt[1, 0].set_yscale('log')
    axt[1, 0].set_xscale('log')
    axt[1, 0].set_ylabel(r'$\pi_k$')
    axt[1, 0].set_xlabel(r'$\lambda_k$ [ns$^{-1}$]')
    axt[1, 0].set_xlim(rmin, rmax)
    axt[1, 0].set_ylim(wmin, wmax)

    axp[1, 0].set_yscale('log')
    axp[1, 0].set_xscale('log')
    axp[1, 0].set_ylabel(r'$\pi_k$')
    axp[1, 0].set_xlabel(r'$\lambda_k$ [ns$^{-1}$]')
    axp[1, 0].set_xlim(rmin, rmax)
    axp[1, 0].set_ylim(wmin, wmax)

    # for individual plots
    ax1a.set_yscale('log')
    ax1a.set_xscale('log')
    ax1a.set_ylabel(r'$\pi_k$')
    ax1a.set_xlabel(r'$\lambda_k$ [ns$^{-1}$]')
    ax1a.set_xlim(rmin, rmax)
    ax1a.set_ylim(wmin, wmax)

    ax1t.set_yscale('log')
    ax1t.set_xscale('log')
    ax1t.set_ylabel(r'$\pi_k$')
    ax1t.set_xlabel(r'$\lambda_k$ [ns$^{-1}$]')
    ax1t.set_xlim(rmin, rmax)
    ax1t.set_ylim(wmin, wmax)

    ax1p.set_yscale('log')
    ax1p.set_xscale('log')
    ax1p.set_ylabel(r'$\pi_k$')
    ax1p.set_xlabel(r'$\lambda_k$ [ns$^{-1}$]')
    ax1p.set_xlim(rmin, rmax)
    ax1p.set_ylim(wmin, wmax)

    for suffix in ['pdf', 'png']:
        basename = (f'basicrta-{gibbs.cutoff}/{gibbs.residue}/weight_vs_rate_'
                    f'results')
        if remove_noise:
            fig1a.savefig(f"{basename}_all_noiserm.{suffix}",
                          bbox_inches='tight')
            fig1t.savefig(f"{basename}_train_noiserm.{suffix}",
                          bbox_inches='tight')
            fig1p.savefig(f"{basename}_validate_noiserm.{suffix}",
                          bbox_inches='tight')
        else:
            fig1a.savefig(f"{basename}_all.{suffix}",
                          bbox_inches='tight')
            fig1t.savefig(f"{basename}_train.{suffix}",
                          bbox_inches='tight')
            fig1p.savefig(f"{basename}_validate.{suffix}",
                          bbox_inches='tight')

    for fig in [fig1a, fig1t, fig1p]:
        plt.close(fig=fig)

    # finish legend for combined plots
    ahandles, aplot_labels = axa[0, 0].get_legend_handles_labels()
    thandles, tplot_labels = axt[0, 0].get_legend_handles_labels()
    phandles, pplot_labels = axp[0, 0].get_legend_handles_labels()

    phs = [plt.plot([], marker="", ls="")[0]]
    pa, pt, pp, pha, pht, php = [], [], [], [], [], []
    for k in range(len(uniq_labels)):
        [pa.append(labl) for labl in aplot_labels[2 * k:2 * (k + 1)]]
        [pt.append(labl) for labl in tplot_labels[k:k + 1]]
        [pp.append(labl) for labl in pplot_labels[k:k + 1]]
        pa.append(f'{len(uniq_labels) - k - 1}')
        pt.append(f'{len(uniq_labels) - k - 1}')
        pp.append(f'{len(uniq_labels) - k - 1}')
        [pha.append(hand) for hand in ahandles[2 * k:2 * k + 2]]
        [pht.append(hand) for hand in thandles[k:k + 1]]
        [php.append(hand) for hand in phandles[k:k + 1]]
        pha.append(phs[0])
        pht.append(phs[0])
        php.append(phs[0])

    ahandles = phs * 3 + pha[::-1]
    thandles = phs * 2 + pht[::-1]
    phandles = phs * 2 + php[::-1]
    aplot_labels = (['Cluster', 'Training', 'Validation'] + pa[::-1])
    tplot_labels = (['Cluster', 'Training'] + pt[::-1])
    pplot_labels = (['Cluster', 'Validation'] + pp[::-1])

    la = figa.legend(ahandles, aplot_labels, loc='lower center',
                     ncols=len(uniq_labels) + 1)
    lt = figt.legend(thandles, tplot_labels, loc='lower center',
                     ncols=len(uniq_labels) + 1)
    lp = figp.legend(phandles, pplot_labels, loc='lower center',
                     ncols=len(uniq_labels) + 1)

    for vpack in la._legend_handle_box.get_children()[1:]:
        for hpack in vpack.get_children()[:1]:
            hpack.get_children()[0].set_width(0)
    for vpack in lt._legend_handle_box.get_children()[1:]:
        for hpack in vpack.get_children()[:1]:
            hpack.get_children()[0].set_width(0)
    for vpack in lp._legend_handle_box.get_children()[1:]:
        for hpack in vpack.get_children()[:1]:
            hpack.get_children()[0].set_width(0)

    figa.tight_layout(rect=(0, 0.1, 1, 1))
    figt.tight_layout(rect=(0, 0.1, 1, 1))
    figp.tight_layout(rect=(0, 0.1, 1, 1))

    figa.subplots_adjust(wspace=0.3, hspace=0.3)
    figt.subplots_adjust(wspace=0.3, hspace=0.3)
    figp.subplots_adjust(wspace=0.3, hspace=0.3)
    for suffix in ['pdf', 'png']:
        basename = f'basicrta-{gibbs.cutoff}/{gibbs.residue}/combined_results'
        if remove_noise:
            figa.savefig(f"{basename}_all_noiserm.{suffix}",
                         bbox_inches='tight')
            figt.savefig(f"{basename}_train_noiserm.{suffix}",
                         bbox_inches='tight')
            figp.savefig(f"{basename}_validate_noiserm.{suffix}",
                         bbox_inches='tight')
        else:
            figa.savefig(f"{basename}_all.{suffix}",
                         bbox_inches='tight')
            figt.savefig(f"{basename}_train.{suffix}",
                         bbox_inches='tight')
            figp.savefig(f"{basename}_validate.{suffix}",
                         bbox_inches='tight')

    if show:
        figa.show()
        figt.show()
        figp.show()
    return all_labels, presorts


def get_code(resname):
    if resname == '-':
        result = '-'
    elif resname == 'X':
        result = 'HSD'
    else:
        result = mda.lib.util.convert_aa_code(resname)
    return result


def get_diffcode(sel, index):
    resid = sel.residues.resids[index]
    resname = sel.residues.resnames[index]
    resletter = mda.lib.util.convert_aa_code(resname)
    return resletter+str(resid)


def get_indices(resnames, sequence):
    indices = []
    inames = 0
    iseq = 0
    while len(indices) < len(sequence):
        if resnames[inames] == sequence[iseq]:
            indices.append(inames)
            iseq += 1
        inames += 1
    return np.asarray(indices)


def get_fa_sel(aln, protA, protB):
    with open(aln) as F:
        names = []
        resids = []
        seqs = []
        tmpseqs = []
        seq = 0
        for line in F:
            if line[0] == '>':
                names.append(line.split('|')[0][1:])
                ress = line.split('/')[1].split('-')
                resids.append([int(i) for i in ress])
                seq += 1
            else:
                tmpseqs.append([seq, line.split('\n')[0]])

    tmpseqs = np.asarray(tmpseqs)
    [seqs.append(''.join(tmpseqs[tmpseqs[:, 0] == k][:, 1])) for k in
     np.unique(tmpseqs[:, 0])]
    seqs = np.asarray(seqs)

    seqA = np.array([i for i in seqs[0]])
    seqB = np.array([i for i in seqs[1]])

    inds = np.where((seqA != '-') & (seqB != '-'))
    selA_mat = protA.residues[inds]
    selB_mat = protB.residues[inds]
    return selA_mat, selB_mat

def get_fa_sel_match(aln, protA, protB):
    with open(aln) as F:
        names = []
        resids = []
        seqs = []
        tmpseqs = []
        seq = 0
        for line in F:
            if line[0] == '>':
                names.append(line.split('|')[0][1:])
                ress = line.split('/')[1].split('-')
                resids.append([int(i) for i in ress])
                seq += 1
            else:
                tmpseqs.append([seq, line.split('\n')[0]])

    tmpseqs = np.asarray(tmpseqs)
    [seqs.append(''.join(tmpseqs[tmpseqs[:, 0] == k][:, 1])) for k in
     np.unique(tmpseqs[:, 0])]
    seqs = np.asarray(seqs)

    seqA = np.array([i for i in seqs[0]])
    seqB = np.array([i for i in seqs[1]])
    match_inds = np.where(seqA == seqB)[0]

    selA_mat = protA.residues[match_inds]
    selB_mat = protB.residues[match_inds]
    return selA_mat, selB_mat

def align_homologues(Areduced, Breduced, aln):
    from MDAnalysis.analysis import align
    uA = mda.Universe(Areduced)
    uB = mda.Universe(Breduced)
                                          
    protA = uA.select_atoms('protein and name CA BB')
    protB = uB.select_atoms('protein and name CA BB')

    selA_mat, selB_mat = get_fa_sel(aln, protA, protB)
    align.alignto(selA_mat, selB_mat)

    uA.atoms.write('Aaligned.pdb')
    uB.atoms.write('Baligned.pdb')

def get_delta_tau(aln, protA, protB, tausA, tausB):
    from MDAnalysis.analysis import align

    residsA = protA.residues.resids
    residsB = protB.residues.resids
    aln_sel = align.fasta2select(aln, is_aligned=True, target_resids=residsB, 
                                 ref_resids=residsA)
    
    selA = protA.select_atoms(aln_sel['reference'])
    selB = protB.select_atoms(aln_sel['mobile'])
    residsA = selA.residues.resids
    residsB = selB.residues.resids

    matchids = np.stack((residsA, residsB)).T
    match_vals = np.array([[tausA[:, 1][tausA[:, 0] == iA][0], 
                            tausB[:, 1][tausB[:, 0] == iB][0], iA, iB] 
                           for iA, iB in matchids if iA in tausA[:, 0] 
                           if iB in tausB[:, 0]])

    delta_tau = -np.diff(match_vals[:, :2]).reshape(len(match_vals),)
    return match_vals[:,2].astype(int), match_vals[:, 3].astype(int), delta_tau

def plot_delta_tau(A, B, dtau, protA, protB, factor=2):
    scale = 1
    rmsd = np.sqrt(np.mean(dtau**2))
    fig, ax = plt.subplots(1, figsize=(4*scale, 3*scale))
    ax.plot(A[dtau > 0], dtau[dtau > 0], '.', color='C0')
    ax.plot(A[dtau < 0], dtau[dtau < 0], '.', color='C3')
    for i, tau in enumerate(dtau):
        if tau >= factor*rmsd:
            resname = protA.select_atoms(f'resid {A[i]}').resnames[0]
            reslet = mda.lib.util.convert_aa_code(resname)
            ax.text(A[i], tau, f'{reslet}{A[i]}')
        elif (tau<0) & (abs(tau) >= factor*rmsd):
            resname = protB.select_atoms(f'resid {B[i]}').resnames[0]
            reslet = mda.lib.util.convert_aa_code(resname)
            ax.text(A[i], tau, f'{reslet}{B[i]}')
        else:
            continue
    ax.xaxis.set_ticks([])
    ax.set_ylabel(r'$\Delta\tau\, [ns]$')
    plt.gca().spines['top'].set_visible(False)
    plt.gca().spines['bottom'].set_visible(False)
    plt.gca().spines['right'].set_visible(False)
    ax.yaxis.set_minor_locator(MultipleLocator(125))
    ax.yaxis.set_major_locator(MultipleLocator(500))
    plt.tight_layout()
    plt.savefig('delta_tau.pdf', bbox_inches='tight')
    plt.savefig('delta_tau.png', bbox_inches='tight')
    plt.show()
