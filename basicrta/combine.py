#!/usr/bin/env python

"""
Command-line interface for combining contact timeseries from multiple repeat runs.

This module provides functionality to combine contact files from multiple 
trajectory repeats, enabling pooled analysis of binding kinetics.
"""

import os
import argparse

class CombineContacts(object):
    """Class to combine contact timeseries from multiple repeat runs.
    
    This class enables pooling data from multiple trajectory repeats and
    calculating posteriors from all data together, rather than analyzing
    each run separately.
    
    :param contact_files: List of contact pickle files to combine
    :type contact_files: list of str
    :param output_name: Name for the combined output file (default: 'combined_contacts.pkl')
    :type output_name: str, optional
    :param validate_compatibility: Whether to validate that files are compatible (default: True)
    :type validate_compatibility: bool, optional
    """
    
    def __init__(self, contact_files, output_name='combined_contacts.pkl', 
                 validate_compatibility=True):
        self.contact_files = contact_files
        self.output_name = output_name
        self.validate_compatibility = validate_compatibility
        
        if len(contact_files) < 2:
            raise ValueError("At least 2 contact files are required for combining")
            
    def _load_contact_file(self, filename):
        """Load a contact pickle file and return data and metadata."""
        if not os.path.exists(filename):
            raise FileNotFoundError(f"Contact file not found: {filename}")
            
        with open(filename, 'rb') as f:
            contacts = pickle.load(f)
            
        metadata = contacts.dtype.metadata
        return contacts, metadata
        
    def _validate_compatibility(self, metadatas):
        """Validate that contact files are compatible for combining."""
        reference = metadatas[0]
        
        # Check that all files have the same atom groups
        for i, meta in enumerate(metadatas[1:], 1):
            # Compare cutoff
            if meta['cutoff'] != reference['cutoff']:
                raise ValueError(f"Incompatible cutoffs: file 0 has {reference['cutoff']}, "
                               f"file {i} has {meta['cutoff']}")
                               
            # Compare atom group selections by checking if resids match
            ref_ag1_resids = set(reference['ag1'].residues.resids)
            ref_ag2_resids = set(reference['ag2'].residues.resids)
            meta_ag1_resids = set(meta['ag1'].residues.resids)
            meta_ag2_resids = set(meta['ag2'].residues.resids)
            
            if ref_ag1_resids != meta_ag1_resids:
                raise ValueError(f"Incompatible ag1 residues between file 0 and file {i}")
            if ref_ag2_resids != meta_ag2_resids:
                raise ValueError(f"Incompatible ag2 residues between file 0 and file {i}")
                
        # Check timesteps and warn if different
        timesteps = [meta['ts'] for meta in metadatas]
        if not all(abs(ts - timesteps[0]) < 1e-6 for ts in timesteps):
            print("WARNING: Different timesteps detected across runs:")
            for i, (filename, ts) in enumerate(zip(self.contact_files, timesteps)):
                print(f"  File {i} ({filename}): dt = {ts} ns")
            print("This may affect residence time estimates, especially for fast events.")
            
    def run(self):
        """Combine contact files and save the result."""
        print(f"Combining {len(self.contact_files)} contact files...")
        
        all_contacts = []
        all_metadatas = []
        
        # Load all contact files
        for i, filename in enumerate(self.contact_files):
            print(f"Loading file {i+1}/{len(self.contact_files)}: {filename}")
            contacts, metadata = self._load_contact_file(filename)
            all_contacts.append(contacts)
            all_metadatas.append(metadata)
            
        # Validate compatibility if requested
        if self.validate_compatibility:
            print("Validating file compatibility...")
            self._validate_compatibility(all_metadatas)
            
        # Combine contact data
        print("Combining contact data...")
        
        # Calculate total size and create combined array
        total_size = sum(len(contacts) for contacts in all_contacts)
        reference_metadata = all_metadatas[0].copy()
        
        # Extend metadata to include trajectory source information
        reference_metadata['source_files'] = self.contact_files
        reference_metadata['n_trajectories'] = len(self.contact_files)
        
        # Determine number of columns (5 for raw contacts, 4 for processed)
        n_cols = all_contacts[0].shape[1]
        
        # Create dtype with extended metadata
        combined_dtype = np.dtype(np.float64, metadata=reference_metadata)
        
        # Add trajectory source column (will be last column)
        combined_contacts = np.zeros((total_size, n_cols + 1), dtype=np.float64)
        
        # Combine data and add trajectory source information
        offset = 0
        for traj_idx, contacts in enumerate(all_contacts):
            n_contacts = len(contacts)
            # Copy original contact data
            combined_contacts[offset:offset+n_contacts, :n_cols] = contacts[:]
            # Add trajectory source index
            combined_contacts[offset:offset+n_contacts, n_cols] = traj_idx
            offset += n_contacts
            
        # Create final memmap with proper dtype
        final_contacts = combined_contacts.view(combined_dtype)
        
        # Save combined contacts
        print(f"Saving combined contacts to {self.output_name}...")
        final_contacts.dump(self.output_name, protocol=5)
        
        print(f"Successfully combined {len(self.contact_files)} files into {self.output_name}")
        print(f"Total contacts: {total_size}")
        print(f"Added trajectory source column (index {n_cols}) for kinetic clustering support")
        
        return self.output_name

def get_parser():
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
    return parser

def main():
    parser = get_parser()
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

if __name__ == "__main__":
    exit(main())    
