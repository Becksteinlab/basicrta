"""
Command line functionality of basicrta.

The `main()` function of this module gets the argument parser from each of the
scripts below and executes the `main()` function of the module called. The
function also collects help from the subparsers and provides it at the command 
line.

Modules callable from the cli: contacts.py, gibbs.py, cluster.py, kinetics.py,
combine.py.
"""

from importlib.metadata import version
import basicrta
import argparse
import subprocess
import importlib
import sys

__version__ = version("basicrta")

# define which scripts can be ran from cli
# can easily add functionality to cli as modules are added
commands = ['contacts', 'gibbs', 'cluster', 'combine', 'kinetics']

def main():
    """ This module provides the functionality for a command line interface for
    basicrta scripts. The scripts available to the cli are:
    
    * contacts.py
    * gibbs.py
    * cluster.py
    * combine.py
    * kinetics.py

    Each script is called and ran using the `main()` function of each module and 
    the parser is passed to the cli using the `get_parser()` function. Any
    module added to the cli needs to have both functions.
    """
    parser = argparse.ArgumentParser(prog='basicrta', add_help=True)
    subparsers = parser.add_subparsers(help="""step in the basicrta workflow to
                                       execute""")
    
    # collect parser from each script in `commands`
    for command in commands:
        subparser = importlib.import_module(f"basicrta.{command}").get_parser()
        subparsers.add_parser(f'{command}', parents=[subparser], add_help=True,
                              description=subparser.description, 
                              conflict_handler='resolve',
                              formatter_class=argparse.ArgumentDefaultsHelpFormatter,
                              help=subparser.description) 
    
    # print subparser help if no arguments given
    if len(sys.argv) == 2 and sys.argv[1] in commands:
        subparsers.choices[f'{sys.argv[1]}'].print_help()
        sys.exit()
    
    # print basicrta help if no subcommand given
    parser.parse_args(args=None if sys.argv[1:] else ['--help'])

    # execute basicrta script
    importlib.import_module(f"basicrta.{sys.argv[1]}").main()

if __name__ == "__main__":
    main()
