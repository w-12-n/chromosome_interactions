# Gene Expression Analysis
This code analyzes DNA data collected by Rao et al. 
The objective is to analyze the spatial organization of DNA within a cell, in order to study how gene expression differs between different cell-types.
We consider the interaction frequencies between all pairs of chromosomes.

We modeled log(1 + interaction frequency) as a Gaussian random variable for each interaction site.
And thus can specify how intensely each pair of chromosomes interacted with each other in the cell.

The results are displayed as a heat map, which specifies the gene expression for a given cell type (in this case, fibrobloast cells).
We may compare heat maps between different cell-types as a gauge of their similarity.

![Heat map for the fibrobloast cells](https://raw.githubusercontent.com/w-12-n/chromosome_interactions/blob/master/heat_map.png)

## Data
The data is split into files, chrX_chrY.txt, which represent sparse matrices of interaction frequencies between chromosome X and chromosome Y in fibrobloast cells. 
The first two columns correspond to the locations on chromosomes X and Y, respectively (in base pairs). 
The third column corresponds to the observed interaction frequency between them. The data's resolution is 250kb.
 
## License
This project is licensed under the MIT License - see the [LICENSE.md](LICENSE.md) file for details
