
#===============================================================================
#===============================================================================
# PROTEIN DEFINITION, CROSS-REFERENCES, ISOFORMS AND DOMAINS
#===============================================================================
#===============================================================================

# Headers and Parameters for Protein object
#===============================================================================

PROTEIN_UNIPROT_DEFINITION_PROPERTY = "PROTEIN_UNIPROT_DEFINITION"

PROTEIN_CLASS = "Protein"

PROTEIN_HEADERS = ["Entry", "Entry name", "Protein names", "Gene names  (primary )", "Gene names  (synonym )", "Organism", "Length", "Fragment", "Cross-reference (Pfam)", "Cross-reference (SMART)"]
PROTEIN_PARAMS = ["Entry", "Entry name", "Protein names", "Gene names  (primary )", "Gene names  (synonym )", "Organism", "Length", "Fragment", "Cross-reference (Pfam)", "Cross-reference (SMART)"]
PROTEIN_COMMENT_CHAR = "#"

# Headers and Parameters for ProteinCrossReference
#===============================================================================

PROTEIN_CROSSREFERENCES_PROPERTY = "PROTEIN_CROSSREFERENCES"

PROTEIN_CROSS_REFERENCE_CLASS = "ProteinCrossReference"

PROTEIN_CROSS_REFERENCE_HEADERS = ["ACC", "DB","X_REFERENCE"]
PROTEIN_CROSS_REFERENCE_PARAMS = ["ACC","DB","X_REFERENCE"]
PROTEIN_CROSS_REFERENCE_COMMENT_CHAR = "#"

# Regular expression and parameters for Isoforms
#===============================================================================

PROTEIN_ISOFORMS_PROPERTY = "PROTEIN_ISOFORMS"

ISOFORM_CLASS = "ProteinIsoform"

ISOFORM_REGULAR_EXPRESSION = ".+\|(?P<ISOFORM_AC>.+)\|.+ Isoform of (?P<PROTEIN_AC>[A-Z0-9]+),.+ OS=(?P<ORGANISM>.+) GN=(?P<GENE_SYMBOL>.+)"
ISOFORM_GROUPS = ["ISOFORM_AC", "PROTEIN_AC"]
ISOFORM_PARAMS = ["ISOFORM_AC", "PROTEIN_AC", "IS_ALTERNATIVE"]
ISOFORM_PARAMS_VALUE_CANONICAL = [None, None, False]
ISOFORM_PARAMS_VALUE_ALTERNATIVE = [None, None, True]
ISOFORM_COMMENT_CHAR = "#"

# Headers and Parameters for ProteinDomain object
#===============================================================================

PROTEIN_DOMAIN_SMART_PROPERTY = "PROTEIN_DOMAIN_SMART"
PROTEIN_DOMAIN_PFAM_PROPERTY = "PROTEIN_DOMAIN_PFAM"

PROTEIN_DOMAIN_CLASS = "ProteinDomain"

PROTEIN_DOMAIN_HEADERS_SMART = ["DOMAIN", "ACC", "DEFINITION", "DESCRIPTION"]
PROTEIN_DOMAIN_PARAM_SMART = ["ACC", "DOMAIN", "DBSOURCE"]
PROTEIN_DOMAIN_DBSOURCE_SMART = "SMART"
PROTEIN_DOMAIN_VALUE_SMART = [None, None, PROTEIN_DOMAIN_DBSOURCE_SMART]

PROTEIN_DOMAIN_HEADERS_PFAM = ["ACC", "DOMAIN", "DEFINITION", "CODE", "DESCRIPTION"]
PROTEIN_DOMAIN_PARAM_PFAM = ["ACC", "DOMAIN", "DBSOURCE"]
PROTEIN_DOMAIN_DBSOURCE_PFAM = "PFAM"
PROTEIN_DOMAIN_VALUE_PFAM = [None, None, PROTEIN_DOMAIN_DBSOURCE_PFAM]

PROTEIN_DOMAIN_COMMENT_CHAR = "-"

#===============================================================================
#===============================================================================
# GENE ONTOLOGY DEFINITION AND ANNOTATIONS
#===============================================================================
#===============================================================================

# Parameters for GeneOntology
#===============================================================================

GENE_ONTOLOGY_DEFINITION_PROPERTY = "GENE_ONTOLOGY_DEFINITION"

GENE_ONTOLOGY_CLASS = "GeneOntology"

GENE_ONTOLOGY_ID_TAG = "id"
GENE_ONTOLOGY_NAME_TAG = "name"
GENE_ONTOLOGY_NAMESPACE_TAG = "namespace"

# Headers and Parameters for ProteinGOAnnotation
#===============================================================================

GENE_ONTOLOGY_ANNOTATION_PROPERTY = "GENE_ONTOLOGY_ANNOTATION"

PROTEIN_GO_ANNOTATION_CLASS = "ProteinGOAnnotation"

PROTEIN_GO_ANNOTATION_HEADERS = ["DB","DB Object ID","DB Object Symbol","Qualifier","GO ID","DB:Reference","Evidence Code","With (or) From","Aspect","DB Object Name","DB Object Synonym","DB Object Type","Taxon(|taxon)","Date","Assigned By","Annotation Extension","Gene Product Form ID"] 
PROTEIN_GO_ANNOTATION_PARAMS = ["DB", "DB Object ID", "GO ID", "Evidence Code"]
PROTEIN_GO_ANNOTATION_COMMENT_CHAR = "!"

#===============================================================================
#===============================================================================
# KEGG_PATHWAY DEFINITION AND ANNOTATIONS
#===============================================================================
#===============================================================================

# Headers and Parameters for KEGGPathway
#===============================================================================

KEGG_PATHWAY_DEFINITION_PROPERTY = "KEGG_PATHWAY_DEFINITION"

KEGG_PATHWAY_CLASS = "KEGGPathway"

KEGG_PATHWAY_HEADERS = ["ID", "Name"]
KEGG_PATHWAY_PARAMS = ["ID","Name"]
KEGG_PATHWAY_COMMENT_CHAR = "#"

# Headers and Parameters for ProteinKEGGannotation
#===============================================================================

KEGG_PATHWAY_ANNOTATION_PROPERTY = "KEGG_PATHWAY_ANNOTATION"

KEGG_PATHWAY_ANNOTATION_CLASS = "ProteinKEGGAnnotation"

KEGG_PATHWAY_ANNOTATION_HEADERS = ["ProteinCrossRef", "KeggID"]
KEGG_PATHWAY_ANNOTATION_PARAMS = ["ProteinCrossRef","KeggID"]
KEGG_PATHWAY_ANNOTATION_COMMENT_CHAR = "#"

#===============================================================================
#===============================================================================
# REACTOME PATHWAY DEFINITION AND ANNOTATIONS
#===============================================================================
#===============================================================================

# Headers and Parameters for ReactomePathway
#===============================================================================

REACTOME_PATHWAY_DEFINITION_PROPERTY = "REACTOME_PATHWAY_DEFINITION"

REACTOME_PATHWAY_CLASS = "ReactomePathway"

REACTOME_PATHWAY_HEADERS = ["ID", "Name", "Organism"]
REACTOME_PATHWAY_PARAMS = ["ID","Name"]
REACTOME_PATHWAY_COMMENT_CHAR = "#"

# Headers and Parameters for ProteinReactomeannotation
#===============================================================================

REACTOME_PATHWAY_ANNOTATION_PROPERTY = "REACTOME_PATHWAY_ANNOTATION"

REACTOME_PATHWAY_ANNOTATION_CLASS = "ProteinReactomeAnnotation"

REACTOME_PATHWAY_ANNOTATION_HEADERS = ["UniprotAC", "ReactomeID", "URL", "ReactomeName", "Tag", "Organism"]
REACTOME_PATHWAY_ANNOTATION_PARAMS = ["UniprotAC", "ReactomeID"]
REACTOME_PATHWAY_ANNOTATION_COMMENT_CHAR = "#"

#===============================================================================
#===============================================================================
# BIOPLEX DEFINITION AND ANNOTATIONS
#===============================================================================
#===============================================================================

# Headers and Parameters for BioplexCluster
#===============================================================================

BIOPLEX_CLUSTER_DEFINITION_PROPERTY = "BIOPLEX_CLUSTER_DEFINITION"

BIOPLEX_CLUSTER_CLASS = "BioplexCluster"

BIOPLEX_CLUSTER_HEADERS = ["ID"]
BIOPLEX_CLUSTER_PARAMS = ["ID"]
BIOPLEX_CLUSTER_COMMENT_CHAR = "#"

# Headers and Parameters for ProteinBioplexAnnotaiton
#===============================================================================

BIOPLEX_ANNOTATION_PROPERTY = "BIOPLEX_ANNOTATION"

BIOPLEX_ANNOTATION_CLASS = "ProteinBioplexAnnotation"

BIOPLEX_ANNOTATION_HEADERS = ["ID", "protein_id"]
BIOPLEX_ANNOTATION_PARAMS = ["ID", "protein_id"]
BIOPLEX_ANNOTATION_COMMENT_CHAR = "#"


#===============================================================================
#===============================================================================
# WAN DEFINITION AND ANNOTATIONS
#===============================================================================
#===============================================================================

# Headers and Parameters for WanCluster
#===============================================================================

WAN_CLUSTER_DEFINITION_PROPERTY = "WAN_CLUSTER_DEFINITION"

WAN_CLUSTER_CLASS = "WanCluster"

WAN_CLUSTER_HEADERS = ["ID"]
WAN_CLUSTER_PARAMS = ["ID"]
WAN_CLUSTER_COMMENT_CHAR = "#"

# Headers and Parameters for ProteinWanAnnotaiton
#===============================================================================

WAN_ANNOTATION_PROPERTY = "WAN_ANNOTATION"

WAN_ANNOTATION_CLASS = "ProteinWanAnnotation"

WAN_ANNOTATION_HEADERS = ["ID", "protein_id"]
WAN_ANNOTATION_PARAMS = ["ID", "protein_id"]
WAN_ANNOTATION_COMMENT_CHAR = "#"


#===============================================================================
#===============================================================================
# CORUM DEFINITION AND ANNOTATIONS
#===============================================================================
#===============================================================================

# Headers and Parameters for CorumCluster
#===============================================================================

CORUM_CLUSTER_DEFINITION_PROPERTY = "CORUM_CLUSTER_DEFINITION"

CORUM_CLUSTER_CLASS = "CorumCluster"

CORUM_CLUSTER_HEADERS = ["ID","Name","Method"]
CORUM_CLUSTER_PARAMS = ["ID","Name","Method"]
CORUM_CLUSTER_COMMENT_CHAR = "#"

# Headers and Parameters for ProteinCorumAnnotaiton
#===============================================================================

CORUM_ANNOTATION_PROPERTY = "CORUM_ANNOTATION"

CORUM_ANNOTATION_CLASS = "ProteinCorumAnnotation"

CORUM_ANNOTATION_HEADERS = ["ID", "protein_id"]
CORUM_ANNOTATION_PARAMS = ["ID", "protein_id"]
CORUM_ANNOTATION_COMMENT_CHAR = "#"


#===============================================================================
#===============================================================================
# INTERACTOME DEFINITION, NETWORK, PARTITION AND ANNOTATIONS
#===============================================================================
#===============================================================================

# Headers and Parameters for ProteinInteraction
#===============================================================================

INTERACTOME_DEFINITION_PROPERTY = "INTERACTOME_DEFINITION"

INTERACTOME_CLASS = "ProteinInteraction"

#INTERACTOME_HEADER = ["#ID(s) interactor A", "ID(s) interactor B", "Alt. ID(s) interactor A", "Alt. ID(s) interactor B", "Alias(es) interactor A", "Alias(es) interactor B", "Interaction detection method(s)", "Publication 1st author(s)", "Publication Identifier(s)", "Taxid interactor A", "Taxid interactor B", "Interaction type(s)", "Source database(s)", "Interaction identifier(s)", "Confidence value(s)"]
#INTERACTOME_PARAMS = ["#ID(s) interactor A", "ID(s) interactor B", "Alt. ID(s) interactor A", "Alt. ID(s) interactor B", "Interaction detection method(s)", "Publication Identifier(s)", "Interaction type(s)", "Source database(s)", "Confidence value(s)"]
INTERACTOME_HEADER = ["#ID(s) interactor A", "Alt. ID(s) interactor A", "ID(s) interactor B", "Alt. ID(s) interactor B", "Interaction detection method(s)", "Publication Identifier(s)", "Source database(s)", "Interaction type(s)", "Confidence value(s)"]
INTERACTOME_PARAMS = ["#ID(s) interactor A", "ID(s) interactor B", "Alt. ID(s) interactor A", "Alt. ID(s) interactor B", "Interaction detection method(s)", "Publication Identifier(s)", "Interaction type(s)", "Source database(s)", "Confidence value(s)"]
INTERACTOME_COMMENT_CHAR = "#"
INTERACTOME_ID_UNIPROTKB_REGEX = "([A-Z0-9]{6})"
INTERACTOME_ID_CROSSREF_REGEX_DICT = {"UniProtKB-ID" : "(^[0-9A-Z]+_[A-Z]+$)", "DIP" : "dip:(DIP-[A-Z0-9]+)", "refseq" : "refseq:(NP_[0-9]+)"}
INTERACTOME_ID_CROSSREF_DB_REGEX_DICT = {"UniProtKB-ID" : "VALUE", "DIP" : "VALUE", "refseq" : "VALUE.%"}
INTERACTOME_PUBMED_REGEX = "^.*pubmed:([a-z0-9]+).*$"
INTERACTOME_SCORE_REGEX = "^\D*([0-9\.]+)\D*$"
INTERACTOME_TYPE_REGEX = "psi-mi:(\"MI:[0-9]+\"\(.*\))"
INTERACTOME_METHOD_REGEX = "psi-mi:(\"MI:[0-9]+\"\(.*\))"
INTERACTOME_SOURCEDB_REGEX = "psi-mi:(\"MI:[0-9]+\"\(.*\))"

# Headers and Parameters for Protein Interaction network
#===============================================================================

INTERACTOME_NETWORK_DEFINITION_PROPERTY = "INTERACTOME_NETWORK_DEFINITION"

INTERACTOME_NETWORK_CLASS = "PPINetworkInteraction"

INTERACTOME_NETWORK_HEADER = ["InteractorA", "InteractorB"]
INTERACTOME_NETWORK_PARAMS = ["InteractorA", "InteractorB", "NetworkName"]
INTERACTOME_NETWORK_COMMENT_CHAR = "#"

# Parameter for NetworkModule
#===============================================================================

INTERACTOME_NETWORK_PARTITION_DEFINITION_PROPERTY = "INTERACTOME_NETWORK_PARTITION_DEFINITION"

INTERACTOME_NETWORK_PARTITION_CLASS = "NetworkModule"

INTERACTOME_NETWORK_PARTITION_CLASS_TAG = ">Class"
INTERACTOME_NETWORK_PARTITION_COMMENT_CHAR = "#"

# Parameter for NetworkModuleAnnotation
#===============================================================================

INTERACTOME_NETWORK_PARTITION_ANNOTATION_PROPERTY = "INTERACTOME_NETWORK_PARTITION_ANNOTATION"

INTERACTOME_NETWORK_PARTITION_ANNOTATION_CLASS = "NetworkModuleAnnotation"

INTERACTOME_NETWORK_PARTITION_ANNOTATION_CLASS_TAG = "[CLASS"
INTERACTOME_NETWORK_PARTITION_ANNOTATION_CLASS_REGEX = "\[CLASS: ([0-9]+)\]"
INTERACTOME_NETWORK_PARTITION_ANNOTATION_PROTEIN_TAG = "PN"
INTERACTOME_NETWORK_PARTITION_ANNOTATION_ANNOTATION_TAG = "CA"
INTERACTOME_NETWORK_PARTITION_ANNOTATION_COMMENT_CHAR = "#"

# Parameter for ProteinCrossReference from protein redundancy
#===============================================================================
INTERACTOME_NETWORK_REDUNDANCY_DEFINITION_PROPERTY = "INTERACTOME_NETWORK_REDUNDANCY_DEFINITION"

INTERACTOME_NETWORK_REDUNDANCY_DEFINITION_CLASS = "ProteinCrossReference"

INTERACTOME_NETWORK_REDUNDANCY_DEFINITION_HEADERS = ["REDUNDANT_AC", "PROTEIN_AC"]
INTERACTOME_NETWORK_REDUNDANCY_DEFINITION_PARAMS = ["PROTEIN_AC","DB_SOURCE","REDUNDANT_AC"]
INTERACTOME_NETWORK_REDUNDANCY_DEFINITION_COMMENT_CHAR = "#"


#===============================================================================
#===============================================================================
# RNA DEFINITION, CROSS-REFERENCES, ISOFORMS AND SUBTYPES
#===============================================================================
#===============================================================================

# Headers and Parameters for RNA object
#===============================================================================

RNA_DEFINITION_PROPERTY = "RNA_DEFINITION"

RNA_CLASS = "RNA"

RNA_HEADERS = ["transcript_ID","parent_gene","peptide_ID","transcript_biotype","transcript_length","transcript_source","transcript_status","transcript_tsl","transcript_gencode_basic","transcript_start","transcript_end","transcript_strand","chromosome_name","percentage_GC_content","description","external_gene_name","external_gene_source","external_transcript_name","external_transcript_source_name"]
RNA_PARAMS = ["transcript_ID","parent_gene","peptide_ID","transcript_biotype","transcript_length","transcript_source","transcript_status","transcript_tsl","transcript_gencode_basic","transcript_start","transcript_end","transcript_strand","chromosome_name","percentage_GC_content","description","external_gene_name","external_gene_source","external_transcript_name","external_transcript_source_name"]
RNA_COMMENT_CHAR = "#"

# RNA types which will reflect the different RNA subtables
RNA_BROAD_TYPES = ["MRNA", "OtherRNA", "LncRNA"]

# MRNA biotypes: all mRNAs are defined as protein_coding. However, not all protein_coding are necessarily mRNAs.
RNA_MRNA_BIOTYPE = ["protein_coding"]
# LncRNA biotypes: Merge from GENCODEv23 and v24 biotypes considered as lncRNAs. See: ftp://ftp.sanger.ac.uk/pub/gencode/Gencode_human/release_23/_README.TXT
RNA_LNCRNA_BIOTYPE = ["3prime_overlapping_ncrna","antisense","bidirectional_promoter_lncrna","known_ncrna","lincRNA","macro_lncRNA","non_coding","processed_transcript","sense_intronic","sense_overlapping","TEC"]
RNA_ALL_BIOTYPE = ["protein_coding","retained_intron","nonsense_mediated_decay","processed_transcript","processed_pseudogene","antisense","unprocessed_pseudogene","lincRNA","unitary_pseudogene","transcribed_unprocessed_pseudogene","transcribed_processed_pseudogene","miRNA","transcribed_unitary_pseudogene","sense_overlapping","pseudogene","sense_intronic","TR_V_gene","3prime_overlapping_ncrna","IG_V_gene","polymorphic_pseudogene","misc_RNA","non_stop_decay","snRNA","IG_V_pseudogene","snoRNA","rRNA","ribozyme","Mt_tRNA","Mt_rRNA","IG_C_gene","IG_J_gene","TR_J_gene","TR_V_pseudogene","TR_C_gene","TR_J_pseudogene","IG_C_pseudogene","IG_D_gene","TR_D_gene","IG_J_pseudogene","TEC","scaRNA","translated_unprocessed_pseudogene","vaultRNA","sRNA","macro_lncRNA"]

RNA_ALL_KW = "allRNAs"
PROT_ALL_KW = "allProts"
PROTEIN_ENSP_XREF_DB = "Ensembl_PRO"
PROTEIN_ENSP_XREF_KW = "ProteinEnsemblPROCrossReference"


# Headers and Parameters for RNACrossReference
#===============================================================================

RNA_CROSS_REFERENCE_PROPERTY = "RNA_CROSS_REFERENCE"

RNA_CROSS_REFERENCE_CLASS = "RNACrossReference"

RNA_CROSS_REFERENCE_HEADERS = ["ACC", "DB","X_REFERENCE"]
RNA_CROSS_REFERENCE_PARAMS = ["ACC","DB","X_REFERENCE"]
RNA_CROSS_REFERENCE_COMMENT_CHAR = "#"

#===============================================================================
#===============================================================================
# PROTEIN RNA INTERACTION (PRI) DEFINITION
#===============================================================================
#===============================================================================

# Headers and Parameters for ProteinRNAInteractionCatRAPID
#===============================================================================

PROTEIN_RNA_INTERACTION_CATRAPID_DEFINITION_PROPERTY = "PROTEIN_RNA_INTERACTION_DEFINITION"

PROTEIN_RNA_INTERACTION_CATRAPID_CLASS = "ProteinRNAInteractionCatRAPID"

PROTEIN_RNA_INTERACTION_CATRAPID_HEADERS = ["interactors", "interactionScore", "otherScore1", "otherScore2"]
PROTEIN_RNA_INTERACTION_CATRAPID_PARAMS = ["interactors", "interactionScore"]
PROTEIN_RNA_INTERACTION_CATRAPID_COMMENT_CHAR = "#"

PROTEIN_RNA_INTERACTION_CATRAPID_MISSING_PROT_KW = "Proteins_not_found"
PROTEIN_RNA_INTERACTION_CATRAPID_MISSING_RNA_KW = "RNAs_not_found"

# Headers and Parameters for InteractingRNA
#===============================================================================

INTERACTING_RNA_DEFINITION_PROPERTY = "INTERACTING_RNA_DEFINITION"
INTERACTING_RNA_DEFINITION_CLASS = "InteractingRNA"

INTERACTING_RNA_DEFINITION_HEADERS = ["ensembl_id"]
INTERACTING_RNA_DEFINITION_PARAMS = ["ensembl_id"]
INTERACTING_RNA_DEFINITION_COMMENT_CHAR = "#"

 
# Headers and Parameters for InteractingProtein
#===============================================================================
 
INTERACTING_PROTEIN_DEFINITION_PROPERTY = "INTERACTING_PROTEIN_DEFINITION"
INTERACTING_PROTEIN_DEFINITION_CLASS = "InteractingProtein"
 
INTERACTING_PROTEIN_DEFINITION_HEADERS = ["uniprotac"]
INTERACTING_PROTEIN_DEFINITION_PARAMS = ["uniprotac"]
INTERACTING_PROTEIN_DEFINITION_COMMENT_CHAR = "#"


#===============================================================================
#===============================================================================
# RNA EXPRESSION DATA
#===============================================================================
#===============================================================================

# Headers and Parameters for RNA Tissue Expression
#===============================================================================

RNA_TISSUE_EXPRESSION_PROPERTY = "RNA_TISSUE_EXPRESSION"

RNA_TISSUE_EXPRESSION_CLASS = "RNATissueExpression"

RNA_TISSUE_EXPRESSION_HEADERS = ["transcript_id", "tissue_name", "rpkm_mean", "rpkm_std", "rpkm_median", "coef_variation", "max"]
RNA_TISSUE_EXPRESSION_PARAMS = ["transcript_id", "tissue_name", "rpkm_mean", "source_db"]
RNA_TISSUE_EXPRESSION_COMMENT_CHAR = "#"

RNA_TISSUE_EXPRESSION_SOURCEDB = "GTEx v6"
RNA_TISSUE_EXPRESSION_VALUE = [None, None, None, RNA_TISSUE_EXPRESSION_SOURCEDB]

