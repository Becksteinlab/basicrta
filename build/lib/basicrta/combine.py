#!/usr/bin/env python

"""
Command-line interface for combining contact timeseries from multiple repeat runs.

This module provides functionality to combine contact files from multiple 
trajectory repeats, enabling pooled analysis of binding kinetics.
"""

import os
import argparse
from basicrta.contacts import CombineContacts


def main():
    """Main function for combining contact files."""
    parser = argparse.ArgumentParser(
        description="Combine contact timeseries from multiple repeat runs. "
                   "This enables pooling data from multiple trajectory repeats "
                   "and calculating posteriors from all data together."
    )
    
    parser.add_argument(
        '--contacts', 
        nargs='+', 
        required=True,
        help="List of contact pickle files to combine (e.g., contacts_7.0.pkl from different runs)"
    )
    
    parser.add_argument(
        '--output', 
        type=str, 
        default='combined_contacts.pkl',
        help="Output filename for combined contacts (default: combined_contacts.pkl)"
    )
    
    parser.add_argument(
        '--no-validate',
        action='store_true',
        help="Skip compatibility validation (use with caution)"
    )
    
    args = parser.parse_args()
    
    # Validate input files exist
    missing_files = []
    for filename in args.contacts:
        if not os.path.exists(filename):
            missing_files.append(filename)
    
    if missing_files:
        print("ERROR: The following contact files were not found:")
        for filename in missing_files:
            print(f"  {filename}")
        return 1
    
    if len(args.contacts) < 2:
        print("ERROR: At least 2 contact files are required for combining")
        return 1
    
    if os.path.exists(args.output):
        print(f"ERROR: Output file {args.output} already exists")
        return 1
    
    try:
        combiner = CombineContacts(
            contact_files=args.contacts,
            output_name=args.output,
            validate_compatibility=not args.no_validate
        )
        
        output_file = combiner.run()
        
        print(f"\nCombination successful!")
        print(f"Combined contact file saved as: {output_file}")
        print(f"\nYou can now use this file with the Gibbs sampler:")
        print(f"  python -m basicrta.gibbs --contacts {output_file} --nproc <N>")
        
        return 0
        
    except Exception as e:
        print(f"ERROR: {e}")
        return 1


if __name__ == '__main__':
    exit(main())
