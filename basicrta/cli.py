"""
basicrta
A package to extract binding kinetics from molecular dynamics simulations
"""

from importlib.metadata import version
from basicrta import *
import argparse
import subprocess
import importlib
import sys

__version__ = version("basicrta")

commands = ['contacts', 'cluster', 'combine', 'kinetics', 'gibbs']

def main():
    parser = argparse.ArgumentParser(prog='basicrta', add_help=True)
    subparsers = parser.add_subparsers()
    
    for command in commands:
        subparser = importlib.import_module(f"basicrta.{command}").get_parser()
        subparsers.add_parser(f'{command}', parents=[subparser], add_help=False)

    args = parser.parse_args(args=None if sys.argv[1:] else ['--help'])

    if len(sys.argv) == 2:
        subparsers.choices[f'{sys.argv[1]}'].print_help()
        sys.exit()
               
    importlib.import_module(f"basicrta.{sys.argv[1]}").main()

if __name__ == "__main__":
    main()
