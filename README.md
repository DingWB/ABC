Forked from [https://github.com/broadinstitute/ABC-Enhancer-Gene-Prediction](https://github.com/broadinstitute/ABC-Enhancer-Gene-Prediction)


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
- EnhancerPredictionsFull_threshold0.025_self_promoter.tsv: take EnhancerPredictionsAllPutative.txt.gz as input and: add non_expressed genes + filter based on ABC score + filter by self promoter 
- EnhancerPredictions_threshold0.025_self_promoter.tsv: subset the columns from EnhancerPredictionsFull_threshold0.025_self_promoter.tsv
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
