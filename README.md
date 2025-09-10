CircleCI [![CircleCI](https://dl.circleci.com/status-badge/img/gh/broadinstitute/ABC-Enhancer-Gene-Prediction.svg?style=svg)](https://app.circleci.com/pipelines/github/broadinstitute/ABC-Enhancer-Gene-Prediction)

> :memo: **Note:** This repository is the version of ABC presented in [Gschwind _et al._ (_BioRxiv_ 2023)](https://doi.org/10.1101/2023.11.09.563812) and is a revamp of the original ABC codebase. To access any previous version of ABC, please visit https://github.com/EngreitzLab/ABC-Enhancer-Gene-Prediction-20250314-archive.
> - For the version of ABC presented in Fulco _et al._ (_Nat. Genet._ 2019)[1], use the [NG2019 branch](https://github.com/EngreitzLab/ABC-Enhancer-Gene-Prediction-20250314-archive/tree/NG2019) of the archived repo
> - For the version of ABC presented in Nasser _et al._ (_Nature_ 2021), use the [master branch](https://github.com/EngreitzLab/ABC-Enhancer-Gene-Prediction-20250314-archive/tree/master) of the archived repo

We have documented the codebase and usage in [Read The Docs](https://abc-enhancer-gene-prediction.readthedocs.io/en/latest/) 

# Activity by Contact Model of Enhancer-Gene Specificity

The Activity-by-Contact (ABC) model predicts which enhancers regulate which genes on a cell type specific basis. This repository contains the code needed to run the ABC model as well as small sample data files, example commands, and some general tips and suggestions.

If you use the ABC model in published research, please cite:

[1] Fulco CP, Nasser J, Jones TR, Munson G, Bergman DT, Subramanian V, Grossman SR, Anyoha R, Doughty BR, Patwardhan TA, Nguyen TH, Kane M, Perez EM, Durand NC, Lareau CA, Stamenova EK, Aiden EL, Lander ES & Engreitz JM. Activity-by-contact model of enhancer–promoter regulation from thousands of CRISPR perturbations. Nat. Genet. 51, 1664–1669 (2019). https://www.nature.com/articles/s41588-019-0538-0

## Contact
Please submit a GitHub issue with any questions or if you experience any issues/bugs.


## Installation
```shell
git clone https://github.com/DingWB/ABC.git
conda env create -y -f ~/Software/ABC-Enhancer-Gene-Prediction/workflow/envs/abcenv.yml --prefix ~/Software/conda/abc
```

## Example
```shell
snakemake -d ~/Software/ABC --snakefile ~/Software/ABC/workflow/run_abc.smk -np
snakemake -d ~/Software/ABC --snakefile ~/Software/ABC/workflow/run_abc.smk -p -j 1
```


## Output
- Neighborhoods/EnhancerList.txt: Candidate enhancer regions with Dnase-seq (or ATAC-seq) and H3K27ac ChIP-seq read counts
- Neighborhoods/GeneList.txt: Dnase-seq (or ATAC-seq) and H3K27ac ChIP-seq read counts on gene bodies and gene promoter regions
- Predictions/EnhancerPredictionsAllPutative.txt.gz: ABC scores for all element-gene pairs. Includes promoter elements and pairs with scores below the threshold. Only includes expressed genes. This file includes both the 'positive' and 'negative' predictions of the model. (use --make_all_putative to generate this file).
- EnhancerPredictionsAllPutativeNonExpressedGenes.txt.gz: Same as above for non-expressed genes. This file is provided for completeness but we generally do not recommend using these predictions.
- Predictions/EnhancerPredictionsFull_threshold0.025_self_promoter.tsv


EnhancerPredictionsFull.txt

| ColumnName                         | example                          | Definition                                                                                                              | Notes                                                         |
| ---------------------------------- | -------------------------------- | ----------------------------------------------------------------------------------------------------------------------- | ------------------------------------------------------------- |
| chr                                | chr1                             | Chromosome of the enhancer and gene                                                                                     |                                                               |
| start                              | 1000000                          | Start coordinate of the enhancer element                                                                                |                                                               |
| end                                | 1000500                          | End coordinate of the enhancer element                                                                                  |                                                               |
| class                              | intergenic                       | Annotation of enhancer element as {promoter,genic,intergenic}                                                           |                                                               |
| name                               | intergenic\|chr1:1000000-1000500 | unique enhancer region identifier                                                                                       |                                                               |
| distance                           | 5000                             | Distance in bp between enhancer and TSS of target gene                                                                  |                                                               |
| isSelfPromoter                     | FALSE                            | Boolean denoting whether element is the promoter of the TargetGene                                                      | Note that EnhancersPredictions.txt does not contain promoters |
| TargetGene                         | TPTEP1                           | Target Gene                                                                                                             |                                                               |
| TargetGeneTSS                      | 1005250                          | Transcription Start Site of Target Gene (used to extract Hi-C Contact)                                                  |                                                               |
| TargetGeneExpression               | 37.96                            | Target Gene Expression                                                                                                  |                                                               |
| TargetGenePromoterActivityQuantile | 0.95                             | Quantile of Activity at Promoter of target gene                                                                         |                                                               |
| hic_contact                        | 0.035                            | K-R normalized Hi-C contacts between element and TargetGene TSS                                                         |                                                               |
| hic_contact_pl_scaled              | 0.035                            | hic_contact scaled by the difference in powerlaw fits between target cell type and reference cell type                  |                                                               |
| hic_pseudocount                    | 0.001                            | pseudocount added to HiC Contact                                                                                        |                                                               |
| hic_contact_pl_scaled_adj          | 0.0351                           | Powerlaw scaled KR Normalized HiC plus pseudocount. This is the Contact used in the ABC Score                           |                                                               |
| activity_base                      | 10.7                             | Geometric mean DHS (or ATAC) and H3K27ac. This is the Activity used in the ABC Score                                    |                                                               |
| ABC.Score.Numerator                | 0.375                            | The numator of the ABC Score. Activity x Contact without dividing by the scores for other elements near the Target Gene |                                                               |

More information: https://github.com/mayasheth/ABC-Enhancer-Gene-Prediction
