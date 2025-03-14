## Genomic_Analysis_of–CARD14_Epidemiological-Variants
### Introduction
Psoriasis is a chronic inflammatory skin disease with a global prevalence of approximately 2–3%, though this rate varies significantly depending on factors such as geographic region, ethnicity, and country, ranging between 0.2% and 10%. Notably, psoriasis is more common among Caucasians, less frequent in Asians, and rarely observed in Africans. Furthermore, studies suggest that psoriasis is more prevalent in high-income regions, particularly Western Europe and North America [1], [2].
The CARD14 protein, primarily expressed in keratinocytes and placental tissues, has been implicated in the development of psoriasis and possibly other skin disorders. Located on chromosome 17, mutations in the CARD14 gene have been recognized as a significant genetic factor contributing to psoriasis susceptibility [3], [4].
This project analyzes genomic data from chromosome 17, obtained from Phase 1 of the 1000 Genomes Project and mapped to the GRCh37 reference genome, to determine whether there is a statistically significant difference in alternate allele frequency between regions with high psoriasis prevalence (EUR – European ancestry, AMR – American ancestry) and regions with low prevalence (AFR – African ancestry, ASN – Asian ancestry) [5].
 
### Methodology [6]
To process the data, two separate Python scripts were developed:
1.	Full gene analysis – Examines the entire genomic region of CARD14 (Chromosome 17: 78,143,811–78,183,130, forward strand).
2.	Exon-specific analysis – Sections the chromosome based on active exon locations (Transcript: ENST00000573882.1).
The CARD14 gene consists of 23 exons, of which 20 are protein-coding exons. This project focuses on analyzing these 20 protein-coding exons to investigate variations in alternate allele frequency.
 
### Observations (Jupyter Notebook)
1. Principal Component Analysis (PCA): No significant clusters were observed in the PCA plot, and no statistical analysis was conducted on this plot.
2. Alternate Allele Frequency Analysis:
 - No statistically significant differences were found in the overall comparison between high- and low-prevalence regions.
 - No statistically significant differences were observed in alternate allele frequency per exon.

[1]	R. Parisi, I. Y. K. Iskandar, E. Kontopantelis, M. Augustin, C. E. M. Griffiths, and D. M. Ashcroft, ‘National, regional, and worldwide epidemiology of psoriasis: systematic analysis and modelling study’, BMJ, vol. 369, p. m1590, May 2020, doi: 10.1136/bmj.m1590.
[2]	E. Christophers, ‘Psoriasis--epidemiology and clinical spectrum’, Clin Exp Dermatol, vol. 26, no. 4, pp. 314–320, Jun. 2001, doi: 10.1046/j.1365-2230.2001.00832.x.
[3]	D. M. Berki et al., ‘Activating CARD14 Mutations Are Associated with Generalized Pustular Psoriasis but Rarely Account for Familial Recurrence in Psoriasis Vulgaris’, Journal of Investigative Dermatology, vol. 135, no. 12, 2015, doi: 10.1038/jid.2015.288.
[4]	M. Mellett et al., ‘CARD14 Gain-of-Function Mutation Alone Is Sufficient to Drive IL-23/IL-17–Mediated Psoriasiform Skin Inflammation In Vivo’, Journal of Investigative Dermatology, vol. 138, no. 9, 2018, doi: 10.1016/j.jid.2018.03.1525.
[5]	‘Human Genome Structural Variation Consortium | 1000 Genomes’. Accessed: Mar. 14, 2025. [Online]. Available: https://www.internationalgenome.org/human-genome-structural-variation-consortium/
[6]	‘Transcript: ENST00000573882.1 (CARD14-001) - Summary - Homo_sapiens - GRCh37 Archive browser 113’. Accessed: Mar. 14, 2025. [Online]. Available: https://grch37.ensembl.org/Homo_sapiens/Transcript/Summary?db=core;g=ENSG00000141527;r=17:78143791-78183130;t=ENST00000573882
