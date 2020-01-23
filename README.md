# MutARNA: mutational analysis and visualization for short and long span interactions of RNAs


## Usage


```
usage: MutaRNA-plot.py [-h] --fasta-wildtype FASTA_WILDTYPE --SNP-tag SNP_TAG
                       [--out-dir OUT_DIR] [--no-global-fold]
                       [--no-local-fold] [--local-W LOCAL_W]
                       [--local-L LOCAL_L] [--global-maxL GLOBAL_MAXL]
                       [--no-SNP-score]

MutaRNA-plot predict and plot local and global base-pair probabilities of
wildtype and mutant RNAs Sample call: "python bin/MutaRNA-plot.py --fasta-
wildtype data/sample0.fa --SNP-tag G3C --out-dir tmp --no-global-fold"

optional arguments:
  -h, --help            show this help message and exit
  --fasta-wildtype FASTA_WILDTYPE
                        Input sequence wildtype in fasta format
  --SNP-tag SNP_TAG     SNP tag e.g. "C3G" for mutation at position 3 from C
                        to G
  --out-dir OUT_DIR     output directory
  --no-global-fold      Do not run (semi-)global fold (semi: max-window
                        1000nt)
  --no-local-fold       Do not run local fold
  --local-W LOCAL_W     Window length for local fold
  --local-L LOCAL_L     Max base-pair interaction span for local fold
  --global-maxL GLOBAL_MAXL
                        Maximum interaction span of global length.
  --no-SNP-score        Do not run SNP structure abberation scores with RNAsnp
                        and remuRNA

```


### Example call

`python bin/MutaRNA-plot.py --fasta-wildtype data/sample0.fa --SNP-tag G3C --out-dir ./results/`

### Output
Outputs are stored under the current directory (by default) or the specified path via `--out-dir` option. 

* Local 

The results under the `local/` directory contain the visualization of base-pairing probabilities as computed by RNAplfold under different conditions. The plots are generated in two file formats (`.png` and `svg`). 

    - `RNA-WILD-circos`: The base pair probabilities of the wild type in circular Circos form. 
    - `RNA-MUTANT-circos`: The base pair probabilities of the mutant in circular Circos form. 
    - `RNA-MUTANT-removed-circos`: The weakened base-pairing potentials in the form of probability difference between WT and mutant.
    - `RNA-MUTANT-introduced-circos`: The increased base-pairings potentials in the form of probability difference between mutant and WT.
    - `RNA-WT-MUT-dotplot`: The base pair probabilities of the wild type (upper right) and mutant (lower left) in dotplot-style matrix format heatmap plots.
    - `RNA-REMOVED-INTRODUCED-dotplot`: The weakened base-pairing potentials in the form of probability difference between WT and mutant (upper right). The increased base-pairings potentials in the form of probability difference between mutant and WT (lower left).
    - `RNA_lunp-ECG` The position-wise unpaired probability, also know as accessibility, of wild type, mutant and their difference.

* Global 

The results under the `global/` directory contain the visualization of base-pairing probabilities as computed by RNAplfold. Similar to the local mode but allowing for large base-pair span over the sequence length in a single window. under different conditions. The plots are generated in two file formats (`.png` and `svg`). 

* Predicted structural impact of the mutation:

    - `remuRNA.csv` : prediction scores by remuRNA. remuRNA score `(H(WT||MUT))`is the relative entropy between the ensemble of structures in wild type versus mutant RNA. The score reflects the changes in the global structure of the RNA.
    - `RNAsnp.csv` : prediction scores by RNAsnp. RNAsnp scores are generated in two modes (`-m 1`, `-m 2`) of  (semi-)global and local folding and based on the two metrics of base-pairing distance (`d_max`) in both modes and correlation coefficient (`r_min`). RNAsnp further computes the significance of score in term of  `p-value`s against a pre-computed table of sequence with similar features.
