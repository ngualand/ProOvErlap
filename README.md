# ProOvErlap - Assessing feature proximity/overlap and testing statistical significance from genomic intervals

Feature overlap is a critical concept in bioinformatics and occurs when two genomic intervals, usually represented as BED files, are located in the same genomic regions. Instead, feature proximity refers to the spatial proximity of genomic elements. For example, promoters typically overlap or are close to the genes they regulate. Overlap and proximity are also important in epigenetic studies. Here, the overlap of regions enriched for specific epigenetic modifications or accessible chromatin can elucidate complex molecular phenotypes. Consequently, the ability to analyze and interpret feature overlap and proximity is essential for understanding the biological processes that contribute to a given phenotype. To address this need, we present a computational method capable of analyzing data represented in the BED format. This method aims to quantitatively assess the degree of proximity or overlap between genomic features and to determine the statistical significance of these events in the context of a non-parametric randomization test. The aim is to understand whether the observed state differs from what would be expected by chance. The method is designed to be easy to use, requiring only a single command line to run, allowing straightforward overlap and proximity analysis. It also provides clear visualizations and publication-quality figures. In conclusion, this study highlights the importance of feature overlap and proximity in epigenetic studies and presents a method to improve the systematic assessment and interpretation of these features. A new resource for identifying biologically significant interactions between genomic features in both healthy and disease states.



Usage:

python3 prooverlap.py --help


Options:

