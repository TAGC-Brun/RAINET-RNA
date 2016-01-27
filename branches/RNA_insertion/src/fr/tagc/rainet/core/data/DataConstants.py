
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
INTERACTOME_ID_UNIPROTKB_REGEX = "(^[A-Z0-9]{6}$)"
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

RNA_HEADERS = ["rna_ID","parent_gene","peptide_ID","transcript_biotype","transcript_length","transcript_source","transcript_status","transcript_tsl","transcript_gencode_basic","transcript_start","transcript_end","transcript_strand","chromosome_name","percentage_GC_content"]
RNA_PARAMS = ["rna_ID","parent_gene","peptide_ID","transcript_biotype","transcript_length","transcript_source","transcript_status","transcript_tsl","transcript_gencode_basic","transcript_start","transcript_end","transcript_strand","chromosome_name","percentage_GC_content"]
RNA_COMMENT_CHAR = "#"

