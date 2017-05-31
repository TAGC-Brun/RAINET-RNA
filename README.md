# RAINET project

RAINET: Towards the understanding of **R**N**A**-protein **I**nteraction **NET**works.

This github page contains the source code associated to the Ribeiro _et al._ publication. All resulting data and plots will be published elsewhere.

This project aimed at identifying human long non-coding RNAs likely to scaffold protein complexes. 


## Authors:

Diogo M. Ribeiro 1, Andreas Zanzoni 1, Lionel Spinelli 1, Riccardo Delli Ponti 2,3, Gian Gaetano Tartaglia 2,3,4, Christine Brun 1,5

## Affiliations:
1 - Aix-Marseille Université, Inserm, TAGC UMR_S1090, Marseille, France

2 - Centre for Genomic Regulation (CRG), The Barcelona Institute of Science and Technology, Dr Aiguader 88, 08003 Barcelona, Spain

3 - Universitat Pompeu Fabra (UPF), 08003 Barcelona, Spain

4 - Institucio Catalana de Recerca i Estudis Avançats (ICREA), 23 Passeig Lluıs Companys, 08010 Barcelona, Spain

5 - CNRS, Marseille, France


## Description:

The main scripts for this project were coded in an object-oriented environment, creating and using a SQLite database to contain and associate information. Analysis leading to publication involved many separate analysis scripts, the main ones will be described here.

### Folder description:
* The 'src' directory contains our Python and R scripts. (main analysis scripts are on: /tagc-rainet-RNA/src/fr/tagc/rainet/core/execution)
* The 'test' directory contains Python unittests of main scripts.
* The 'resources' directory contains paths to input files used for main scripts.
* The 'documentation' directory contains the doxygen script documentation and the RAINET database schema.

### Pipeline description:

* RAINET DB preparation:
  * build_data_freeze.sh (get protein-related data)
  * rna_biomart.R (get rna-related data)
  * ProcessGTExData.py (process expression data)
  * Macromolecular complexes (process protein complex data):
    * parse_bioplex_table.py 
    * parse_corum_table.py
    * parse_wan_table.py
  * ReadCatrapid.py (process catRAPID software data (licenced))
  * InsertionStrategy.py (for creating a RAINET database)
  * AnalysisStrategy.py (for filtering / stats of RAINET database)
* Enrichment analysis
  * Main scripts:
    * EnrichmentAnalysisStategy.py (produces main RNA-complex enrichment results)
    * FilterEnrichmentResults.py (for filtering / stats Enrichment analysis results)
  * Enrichment-Disease analysis:
    * ParseLnc2cancer.py (for processing data from lnc2cancer database)
    * ParseLncrnadisease.py (for processing data from lncRNADisease database)
    * OMIM_biomart.R (for retrieving OMIM-related data)
    * OMIMProteinDisease.py (for creating protein-disease correspondence)
    * CommonLncRNAProteinDisease.py (for matching lncRNA and protein disease)
* Other post-analysis:
  * PrioritizeCandidates.py (for selecting enrichments with known interactions)
  * LncRNAGroupOddsRatio.py (for evaluating overlap of groups of lncRNAs against functional lncRNAs)
  * ReadGulko.java and gulko2015_fitcons.R (for evaluating fitCons scores)
  * ComplexDatasetOverlap.py (for evaluating overlap between complex datasets)

## Supported OS and required software:

* Linux
* catRAPID omics (licenced software: http://s.tartaglialab.com/page/catrapid_group)
* Bash 4.3
* SQLite 2.8
* R v3.0.2
* Python 2.7
* python packages: sqlalchemy, pandas, numpy, igraph, pickle

## Contacts:
christine-g.brun@inserm.fr

gian.tartaglia@crg.eu


Diogo Ribeiro,

31 May, 2017
