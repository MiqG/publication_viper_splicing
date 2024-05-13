# A recurrent cancer splicing program

## Outline
1. `program_definition.smk`: define cancer program and survival analysis.
2. `immune_evasion.smk`: study how program recapitulates immune evasion.
3. `proliferation.smk`: study how program recapitulates cell proliferation.
4. `tumorigenesis.smk`: study how program recapitulates tumorigenesis.        
        
## Important remarks

Make sure to have all packages required to run the scripts.

## Recommendations
Run the workflow using
```
snakemake --cores 6
```
In case you want to run the rules on your cluster, refer to [snakemake documentation](https://snakemake.readthedocs.io/en/stable/executing/cluster.html).