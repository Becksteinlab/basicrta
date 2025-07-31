"""
basicrta
A package to extract binding kinetics from molecular dynamics simulations
"""

# Add imports here
from importlib.metadata import version
from basicrta import *
import argparse
import subprocess

__version__ = version("basicrta")

commands = ['contacts', 'cluster', 'combine', 'kinetics', 'gibbs']

def main():
    parser = argparse.ArgumentParser(prog='basicrta')
    #parser.add_argument('command', help='Step in workflow to execute', nargs='+')
    subparsers = parser.add_subparsers(dest='command')

    
    parserA = subparsers.add_parser('contacts', help='Ahelp')
    parserA.add_argument('--top', type=str, help='Topology')
    parserA.add_argument('--traj', type=str, help='Trajectory')
    parserA.add_argument('--sel1', type=str, help='First selection (group for'
                         'which tau is to be calculated)')
    parserA.add_argument('--sel2', type=str, help='Second selection (group of'
                         'interest in interactions with first selection)')
    parserA.add_argument('--cutoff', type=float, help='Value to use (in A) for'
                         'the maximum separation distance that constitutes a'
                         'contact.')
    parserA.add_argument('--nproc', type=int, default=1, help='Number of'
                         'processes to use')
    parserA.add_argument('--nslices', type=int, default=100, help='Number of'
                         'trajectory segments to use (if encountering a'
                         'memoryerror, try using a greater value)')

    parserB = subparsers.add_parser('combine', add_help=True, help='Bhelp')
    parserB.add_argument('--contacts', nargs='+', required=True,
                         help="List of contact pickle files to combine (e.g.," 
                         "contacts_7.0.pkl from different runs)")
    parserB.add_argument( '--output', type=str, default='combined_contacts.pkl',
                         help="Output filename for combined contacts (default:"
                         "combined_contacts.pkl)")
    parserB.add_argument( '--no-validate', action='store_true', help="Skip"
                         "compatibility validation (use with caution)")

    parserC = subparsers.add_parser('cluster', help='Chelp')
    parserC.add_argument('--nproc', type=int, default=1)
    parserC.add_argument('--cutoff', type=float)
    parserC.add_argument('--niter', type=int, default=110000)
    parserC.add_argument('--prot', type=str, default=None, nargs='?')
    parserC.add_argument('--label-cutoff', type=float, default=3,
                         dest='label_cutoff',
                         help='Only label residues with tau > '
                         'LABEL-CUTOFF * <tau>. ')
    parserC.add_argument('--structure', type=str, nargs='?')
    # use  for default values
    parserC.add_argument('--gskip', type=int, default=1000, 
                         help='Gibbs skip parameter for decorrelated samples;'
                         'default from https://pubs.acs.org/doi/10.1021/acs.jctc.4c01522')
    parserC.add_argument('--burnin', type=int, default=10000, 
                         help='Burn-in parameter, drop first N samples as equilibration;'
                         'default from https://pubs.acs.org/doi/10.1021/acs.jctc.4c01522')

    parserD = subparsers.add_parser('gibbs', help='Dhelp')
    parserD.add_argument('--contacts')
    parserD.add_argument('--resid', type=int, default=None)
    parserD.add_argument('--nproc', type=int, default=1)
    parserD.add_argument('--niter', type=int, default=110000)
    parserD.add_argument('--ncomp', type=int, default=15)

    parserE = subparsers.add_parser('kinetics', help='Ehelp')
    parserE.add_argument("--gibbs", type=str)
    parserE.add_argument("--contacts", type=str)
    parserE.add_argument("--top_n", type=int, nargs='?', default=None)
    parserE.add_argument("--step", type=int, nargs='?', default=1)
    parserE.add_argument("--wdensity", action='store_true')
    
    args = parser.parse_args()
    keys, values = vars(args).keys(), vars(args).values()
    inarr = [[f"--{key}", f"{value}"] for key, value in zip(keys, values) if key!='command']
    inlist = [aset for alist in inarr for aset in alist]
    subprocess.run(['python',
                    f'/home/r2/opt/basicrta/basicrta/{args.command}.py'] + 
                   inlist)

if __name__ == "__main__":
    main()
