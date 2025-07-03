# Combining Contact Timeseries from Multiple Repeats

This functionality allows you to combine contact timeseries from multiple repeat runs to analyze pooled data and calculate posteriors from all data together, rather than analyzing each run separately.

## Use Cases

- Analyze data from multiple repeat simulations together
- Pool binding events from multiple trajectories for better statistics
- Calculate combined residence time distributions and confidence intervals

## Usage

### 1. Command Line Interface

After generating contact files for each repeat run individually:

```bash
# Run contact analysis for each repeat
python -m basicrta.contacts --top sys.pdb --traj run1.xtc --sel1 "protein" --sel2 "resname CHOL" --cutoff 7.0
mv contacts_7.0.pkl contacts_run1_7.0.pkl

python -m basicrta.contacts --top sys.pdb --traj run2.xtc --sel1 "protein" --sel2 "resname CHOL" --cutoff 7.0  
mv contacts_7.0.pkl contacts_run2_7.0.pkl

python -m basicrta.contacts --top sys.pdb --traj run3.xtc --sel1 "protein" --sel2 "resname CHOL" --cutoff 7.0
mv contacts_7.0.pkl contacts_run3_7.0.pkl

# Combine the contact files
python -m basicrta.combine --contacts contacts_run1_7.0.pkl contacts_run2_7.0.pkl contacts_run3_7.0.pkl --output combined_contacts_7.0.pkl

# Run Gibbs sampler on combined data
python -m basicrta.gibbs --contacts combined_contacts_7.0.pkl --nproc 5
```

### 2. Python API

```python
from basicrta.contacts import CombineContacts

# Combine contact files
combiner = CombineContacts(
    contact_files=['contacts_run1_7.0.pkl', 'contacts_run2_7.0.pkl', 'contacts_run3_7.0.pkl'],
    output_name='combined_contacts_7.0.pkl'
)

output_file = combiner.run()
print(f"Combined contacts saved to: {output_file}")
```

## Features

### Compatibility Validation

The combiner automatically validates that contact files are compatible:

- **Same cutoff distance**: All files must use the same cutoff
- **Same atom groups**: Protein and ligand selections must match
- **Timestep warnings**: Warns if different timestep values are detected across runs

### Metadata Preservation

Combined files preserve and extend metadata:

- Original trajectory information for each source file
- Number of trajectories combined
- Source file tracking for potential kinetic clustering

### Trajectory Source Tracking

Each contact in the combined file includes trajectory source information:
- Original contact data columns preserved
- Additional column with trajectory index for kinetic clustering support

## Limitations

### Kinetic Clustering

Kinetic clustering is **not yet supported** for combined contact data. The code will:

1. Issue warnings when loading combined contact files
2. Raise a clear error if clustering is attempted
3. Suggest alternatives for kinetic clustering analysis

For kinetic clustering, analyze each trajectory separately or implement the extended clustering algorithm that uses trajectory source information.

### Different Trajectory Properties

- **Different timesteps**: The combiner warns about different timestep values but proceeds. This may affect residence time estimates for fast events.
- **Different particle counts**: Unlike trajectory concatenation, this approach handles trajectories with different numbers of particles correctly.

## Error Handling

The combiner includes comprehensive error checking:

```bash
# Missing files
python -m basicrta.combine --contacts file1.pkl missing_file.pkl
# ERROR: Contact file not found: missing_file.pkl

# Incompatible cutoffs  
python -m basicrta.combine --contacts contacts_7.0.pkl contacts_8.0.pkl
# ERROR: Incompatible cutoffs: file 0 has 7.0, file 1 has 8.0

# Skip validation (use with caution)
python -m basicrta.combine --contacts file1.pkl file2.pkl --no-validate
```

## Output Format

Combined contact files:
- Maintain the same format as individual contact files
- Include extended metadata with source tracking
- Add trajectory source column (last column) for each contact
- Can be used directly with existing Gibbs sampler workflow

The Gibbs sampler will process combined files normally but issue warnings about kinetic clustering limitations.