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

commands = ['contacts', 'gibbs', 'cluster', 'combine', 'kinetics']
parser_help = '''
Step in the basicrta workflow to execute. 
'''

def main():
    parser = argparse.ArgumentParser(prog='basicrta', add_help=True)
    subparsers = parser.add_subparsers(help=parser_help)
    
    for command in commands:
        subparser = importlib.import_module(f"basicrta.{command}").get_parser()
        subparsers.add_parser(f'{command}', parents=[subparser], add_help=True,
                              description=subparser.description, 
                              conflict_handler='resolve',
                              formatter_class=argparse.ArgumentDefaultsHelpFormatter,
                              help=subparser.description) 
    
    if len(sys.argv) == 2 and sys.argv[1] in commands:
        subparsers.choices[f'{sys.argv[1]}'].print_help()
        sys.exit()
    
    parser.parse_args(args=None if sys.argv[1:] else ['--help'])
    importlib.import_module(f"basicrta.{sys.argv[1]}").main()

if __name__ == "__main__":
    main()
