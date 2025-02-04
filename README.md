# DISCO

Decomposition Into Single-COpy gene trees ([DISCO](https://doi.org/10.1093/sysbio/syab070)) is a method for decomposing multi-copy gene-family trees while attempting to preserve orthologs and discard paralogs. These single-copy gene trees can be subsequently used by methods that can estimate species trees from single-copy gene trees such as [ASTRAL](https://github.com/smirarab/ASTRAL) or [ASTRID](https://github.com/pranjalv123/ASTRID) in order to obtain an accurate estimation of the species tree. Additionally, DISCO can be paired with concatenation analysis using the script `ca_disco.py`. 

**NOTE:** For species tree estimation default settings are recommended; however, for orthology detection using `-m 2` is recommended so small groups are retrieved. It is also highly recommended that you use the most recent version of DISCO, as it deals with some limitations of TreeSwift.

# DISCO-R

DISCO-R (DISCO with a Refined pruning strategy) is a modification of DISCO that provably retains more speciation-driven quartets than DISCO. This repository contains the variants of DISCO used in the experiments of Hakim et al. 


## Versions

DISCO:
- v1.0 (May 21st 2021): Initial Release
- v1.1 (July 19th, 2021): Added `ca_disco.py` and adjustable filter (instead of just filtering all trees with less than 4 taxa)
- v1.2 (June 8th, 2022): Tree rerooting fix (addresses issue with treeswift rooting with no edge lengths)
- v1.3 (July 8th, 2022): Allows for more flexibility with leaf labels
- v1.3.1 (October 31st, 2022): Fixes issue where extra newline character at the end of the file causes a crash

DISCO-R
- v1.0 (Feb 7th 2024): Readme Updated

**Note:** At present, it is recommended to use the latest version if you want to run the DISCO algorithm as described in Willson et al. 2021, as using earlier versions can cause results inconsistent with the DISCO algorithm in rare cases. 

## Citation

If you use DISCO, please cite:
```
@article{willson2022disco,
  title={DISCO: Species tree inference using multicopy gene family tree decomposition},
  author={Willson, James and Roddur, Mrinmoy Saha and Liu, Baqiao and Zaharias, Paul and Warnow, Tandy},
  journal={Systematic biology},
  volume={71},
  number={3},
  pages={610--629},
  year={2022},
  publisher={Oxford University Press}
}
```

If you use DISCO-R, please cite:
```
@article {Hakim2023.12.05.570122,
	author = {Sheikh Azizul Hakim and MD Rownok Zahan Ratul and Md. Shamsuzzoha Bayzid},
	title = {wQFM-DISCO: DISCO-enabled wQFM improves phylogenomic analyses despite the presence of paralogs},
	elocation-id = {2023.12.05.570122},
	year = {2023},
	doi = {10.1101/2023.12.05.570122},
	publisher = {Cold Spring Harbor Laboratory},
	journal = {bioRxiv}
}
```

## Algorithm

Given a list of multi-copy gene trees, DISCO does the following for each tree:

1. Root the tree and tag each internal vertex as either a duplication event or a speciation event in such a way that minimizes the total number of duplications and losses. We do this with the ASTRAL-Pro rooting and tagging algorithm ([Zhang et. al. 2020](https://doi.org/10.1093/molbev/msaa139)).
2. Decompose gene tree by splitting off the smallest subtree under every vertex tagged as a duplication from the leaves to the root until all duplication events are resolved; it returns the set of single-copy trees produced.

Please check ``variants.pdf`` for the DISCO variants introduced by Hakim et al. 

## Dependencies

- Python 3
- [TreeSwift](https://github.com/niemasd/TreeSwift)

Treeswift can be installed with: `pip install treeswift`

## Usage

### disco.py

**Input**: File containing list of multi-copy trees in newick format

**Output**: File containing resulting list of single-copy trees after decomposition in newick format

```
python3 disco.py -i <input_file> -o <ouput_file> -d <delimiter>
```

#### Arguments

##### Required

```
-i, --input           Input newick tree file
```

##### Optional

```
-o, --output          Output newick tree file
-d, --delimiter       Delimiter separating species name from rest of leaf label
-s, --single-tree     Output only single large tree per gene tree 
-n, --nth-delimiter   Split on Nth delimiter (only works with -d)
-m, --minimum         Minimum number of taxa required for tree to be outputted
-v, --verbose         Enable verbose output
--keep-labels         Keep original leaf labels instead of using species name
--no-decomp           Outputs rooted trees without decomposition
--outgroups           Write outgroups (including ties) to txt file
--remove-in-paralogs  Remove in-paralogs before rooting/scoring
```

#### Example

```cmd
python3 tag_decomp.py -i example/gtrees-mult.trees
```

### ca_disco.py

**Input**: File containing list of multi-copy trees in newick format and set of alignment files in phylip format corresponding to the gene families.

**Output**: Concatenated alignment file in the phylip format

```
python3 ca_disco.py -i <input_trees> -a <alignments_list> -t <taxa_list> -o <output> -d <delimiter> -m <n> 
```

`disco.py` must be present in the same directory as `ca_disco.py` in order for it to run. Also, unlike `disco.py`, it is necessary for the input newick trees given to `ca_disco.py` to have unique leaf labels where the taxon name comes first and is separated from the rest of the name by some delimiter. 

#### Arguments

##### Required

```
-i, --input           Input newick tree file
-a, --alignment       Text file containing paths to alignment files
-t, --taxonset        Text file containing taxa list
-o, --output          Output concatenated alignment file
```

##### Optional
```
-m, --filter          Minimum number of taxa required sequence group to be included
-d, --delimiter       Delimiter separating species name from rest of leaf label
```

#### Example

```cmd
python3 ca_disco.py -i example/g_100.trees -o example.phy -a example/seq_list.txt -t example/taxa_list.txt
```
