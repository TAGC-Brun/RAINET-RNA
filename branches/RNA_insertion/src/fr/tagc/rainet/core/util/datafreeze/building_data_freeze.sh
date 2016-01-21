#!/bin/bash

# Important note : this script is planned to be executed on the Brun server (192.168.1.15) using the rainet user.
# 
# To use this script:
#
# PHASE 0
#  - declare the root folder you want the data copied to in the ROOT_DATA_FOLDER variable below.
#  - check the MoonGO script folder in the ROOT_MOONGO_FOLDER variable 
#  - check the ClassAnnotation scripts folders in the ROOT_CLASS_ANNOTATION_FOLDER variable
#  - check the CD-Hit folder in the CDHIT_FOLDER variable
#  - Copy this script both with the 'update_uniprot.sh' and 'extract_species_from_flat.pl' scripts to the ROOT_DATA_FOLDER
#
# PHASE 1
#  - Launch this script using the command : ./building_datafreeze.sh 2>datafreeze_step1.log
#  - Select the option 1 in the menu
#  - Wait until the end of execution. It could be very long since all the data will be downloaded and all the interactomes builded with MoonGO
# 		Note that:
# 			* A folder for each species (human, fly, yeast) will be created to contain the data in the ROOT_DATA_FOLDER. 
#  			* the files downloaded and built by the 'update_uniprot.sh' and 'extract_species_from_flat.pl' scripts (trembl sequence flat files used by MoonGO)
#    			will remain into the ROOT_DATA_FOLDER and will also be copied to the ROOT_MOONGO_FOLDER folder.
#  			* the MoonGO scripts creating the interactome and PPI files for each species will 
#    			be launched. Their results will be copied back to the correct species data folder in ROOT_DATA_FOLDER.
#
# PHASE 2
#  - For each species, go to the species data folder and open the $SPECIES.binary.nr0.95.gr files in cytoscape.
#  - Remove the part of the graph not connected to the main part
#  - Save the graph in a file (with two column format) called $SPECIES.binary.nr0.95.connected.gr
#
# PHASE 3
#  - Launch again this script using the command : ./building_datafreeze.sh 2>datafreeze_step2.log
#  - Select the option 2 in the menu
#  - Wait until the end. It could be very long since OCG will be launch on all connected PPI networks
#
# PHASE 4
# Once everything is done, the user can verify that the 'insertion_$SPECIES.ini' files in the resources folder have the correct
# file pathes and names. He can then use Rainet2toad to insert the data into the database (one database per species).
#
# python Rainet.py Insertion -s $SPECIES -d /tmp/Rainet2toad/rainet2016.$SPECIES.db.sqlite -i /big-vol/Rainet2toad/resources/insertion_$SPECIES.ini
#

ROOT_DATA_FOLDER=/big-vol/Rainet2toad/data2016

HUMAN_DATA_FOLDER=human
DROSOPHILA_DATA_FOLDER=drome
YEAST_DATA_FOLDER=yeast

ROOT_MOONGO_FOLDER=/big-vol/MoonGO/scripts
ROOT_CLASS_ANNOTATION_FOLDER=/big-vol/ClassAnnotation/
CDHIT_FOLDER=/home/moongo/Software/cd-hit-v4.6.4-2015-0603
 
#############################################################################################
#############################################################################################
#### START OF FIRST PART FUNCTION: download_and_prepare_data()
#############################################################################################
#############################################################################################
download_and_prepare_data(){
	
# ------------------------------------------------------------------------
# HUMAN
# ------------------------------------------------------------------------
	
	echo '#######################################'
	echo '### HUMAN                      ########'
	echo '#######################################'
	
	mkdir -p $ROOT_DATA_FOLDER/$HUMAN_DATA_FOLDER
	cd $ROOT_DATA_FOLDER/$HUMAN_DATA_FOLDER
	
#
# [PROTEIN]
#
	
# PROTEIN_UNIPROT_DEFINITION
	echo '##### Downloading PROTEIN_UNIPROT_DEFINITION'
	wget 'http://www.uniprot.org/uniprot/?query=organism:"Homo sapiens (Human) [9606]"+AND+proteome:up000005640&columns=id,entry name,protein names,genes(PREFERRED),genes(ALTERNATIVE),organism,length,fragment,database(Pfam),database(SMART)&format=tab' -O human_uniprot_protein_list.txt 
	
# PROTEIN_CROSSREFERENCES
	echo '##### Downloading PROTEIN_CROSSREFERENCES'
	wget ftp://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/idmapping/by_organism/HUMAN_9606_idmapping.dat.gz
	
# PROTEIN_ISOFORMS
	echo '##### Downloading PROTEIN_ISOFORMS'
	wget ftp://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/reference_proteomes/Eukaryota/UP000005640_9606_additional.fasta.gz
	
# PROTEIN_DOMAIN_SMART
	echo '##### Downloading PROTEIN_DOMAIN_SMART'
	wget http://smart.embl.de/smart/descriptions.pl -O all_protein_domains_smart.txt
	
# PROTEIN_DOMAIN_PFAM
	echo '##### Downloading PROTEIN_DOMAN PFAM'
	wget ftp://ftp.ebi.ac.uk/pub/databases/Pfam/current_release/Pfam-A.clans.tsv.gz
	
#
# [GENE_ONTONLOGY]
#
# GENE_ONTOLOGY_DEFINITION
	echo '##### Downloading GENE_ONTOLOGY_DEFINITION'
	wget http://purl.obolibrary.org/obo/go/go-basic.obo
	
# GENE_ONTOLOGY_ANNOTATION
	echo '##### Downloading GENE_ONTOLOGY_ANNOTATION'
	wget ftp://ftp.ebi.ac.uk/pub/databases/GO/goa/HUMAN/gene_association.goa_human.gz
	
#
# [KEGG_PATHWAY]
#
	
# KEGG_PATHWAY_DEFINITION
	echo '##### Downloading KEGG_PATHWAY_DEFINITION'
	wget http://rest.kegg.jp/list/pathway/hsa -O human_kegg_pathway_definitions.txt
	
# KEGG_PATHWAY_ANNOTATION
	echo '##### Downloading KEGG_PATHWAY_ANNOTATION'
	wget http://rest.kegg.jp/link/pathway/hsa -O human_kegg_pathway_annotations.txt
	
#
# [REACTOME PATHWAY]
# 
	
# REACTOME_PATHWAY_DEFINITION
	echo '##### Downloading REACTOME_PATHWAY_DEFINITION'
	wget http://www.reactome.org/download/current/ReactomePathways.txt -O all_reactome_pathway_definitions.txt
	
# REACTOME_PATHWAY_ANNOTATION
	echo '##### Downloading REACTOME_PATHWAY_ANNOTATION'
	wget http://www.reactome.org/download/current/UniProt2Reactome.txt -O all_reactome_pathway_annotations.txt
	
#
# [INTERACTOME]
#
	
# See INTERACTOME section below in the "All three species" part
	
#
# Unzip all compressed files
#
	gunzip $ROOT_DATA_FOLDER/$HUMAN_DATA_FOLDER/*.gz
	
# ------------------------------------------------------------------------
# DROSOPHILA
# ------------------------------------------------------------------------
	
	echo '#######################################'
	echo '### DROSOPHILA                 ########'
	echo '#######################################'
	
	mkdir -p $ROOT_DATA_FOLDER/$DROSOPHILA_DATA_FOLDER
	cd $ROOT_DATA_FOLDER/$DROSOPHILA_DATA_FOLDER
	
#
# [PROTEIN]
#
	
# PROTEIN_UNIPROT_DEFINITION
	echo '##### Downloading PROTEIN_UNIPROT_DEFINITION'
	wget 'http://www.uniprot.org/uniprot/?query=organism:"Drosophila melanogaster (Fruit fly) [7227]" AND proteome:up000000803&columns=id,entry name,protein names,genes(PREFERRED),genes(ALTERNATIVE),organism,length,fragment,database(Pfam),database(SMART)&format=tab' -O drome_uniprot_protein_list.txt
	
# PROTEIN_CROSSREFERENCES
	echo '##### Downloading PROTEIN_CROSSREFERENCES'
	wget ftp://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/idmapping/by_organism/DROME_7227_idmapping.dat.gz
	
# PROTEIN_ISOFORMS
	echo '##### Downloading PROTEIN_ISOFORMS'
	wget ftp://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/reference_proteomes/Eukaryota/UP000000803_7227_additional.fasta.gz
	
# PROTEIN_DOMAIN_SMART
	echo '##### Downloading PROTEIN_DOMAIN_SMART'
	wget http://smart.embl.de/smart/descriptions.pl -O all_protein_domains_smart.txt
	
# PROTEIN_DOMAIN_PFAM
	echo '##### Downloading PROTEIN_DOMAIN_PFAM'
	wget ftp://ftp.ebi.ac.uk/pub/databases/Pfam/current_release/Pfam-A.clans.tsv.gz
	
#
# [GENE_ONTONLOGY]
#
	
# GENE_ONTOLOGY_DEFINITION
	echo '##### Downloading GENE_ONTOLOGY_DEFINITION'
	wget http://purl.obolibrary.org/obo/go/go-basic.obo
	
# GENE_ONTOLOGY_ANNOTATION
	echo '##### Downloading GENE_ONTOLOGY_ANNOTATION'
	wget ftp://ftp.ebi.ac.uk/pub/databases/GO/goa/FLY/gene_association.goa_fly.gz
	
#
# [KEGG_PATHWAY]
#
	
# KEGG_PATHWAY_DEFINITION
	echo '##### Downloading KEGG_PATHWAY_DEFINITION'
	wget http://rest.kegg.jp/list/pathway/dme -O drome_kegg_pathway_definitions.txt
	
# KEGG_PATHWAY_ANNOTATION
	echo '##### Downloading KEGG_PATHWAY_ANNOTATION'
	wget http://rest.kegg.jp/link/pathway/dme -O drome_kegg_pathway_annotations.txt
	
#
# [REACTOME PATHWAY]
# 
	
# REACTOME_PATHWAY_DEFINITION
	echo '##### Downloading REACTOME_PATHWAY_DEFINITION'
	wget http://www.reactome.org/download/current/ReactomePathways.txt -O all_reactome_pathway_definitions.txt
	
# REACTOME_PATHWAY_ANNOTATION
	echo '##### Downloading REACTOME_PATHWAY_ANNOTATION'
	wget http://www.reactome.org/download/current/UniProt2Reactome.txt -O all_reactome_pathway_annotations.txt
	
#
# [INTERACTOME]
#
	
# See INTERACTOME section below in the "All three species" part 
	
#
# Unzip all compressed files
#
	gunzip $ROOT_DATA_FOLDER/$DROSOPHILA_DATA_FOLDER/*.gz
	
# ------------------------------------------------------------------------
# YEAST
# ------------------------------------------------------------------------
	
	mkdir -p $ROOT_DATA_FOLDER/$YEAST_DATA_FOLDER
	cd $ROOT_DATA_FOLDER/$YEAST_DATA_FOLDER
	
	
	echo '#######################################'
	echo '### YEAST                      ########'
	echo '#######################################'
	
#
# [PROTEIN]
#
	
# PROTEIN_UNIPROT_DEFINITION
	echo '##### Downloading PROTEIN_UNIPROT_DEFINITION'
	wget 'http://www.uniprot.org/uniprot/?query=organism:"Saccharomyces cerevisiae (strain ATCC 204508 / S288c) (Bakers yeast) [559292]" AND proteome:up000002311&columns=id,entry name,protein names,genes(PREFERRED),genes(ALTERNATIVE),organism,length,fragment,database(Pfam),database(SMART)&format=tab' -O yeast_uniprot_protein_list.txt
	
# PROTEIN_CROSSREFERENCES
	echo '##### Downloading PROTEIN_CROSSREFERENCES'
	wget ftp://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/idmapping/by_organism/YEAST_559292_idmapping.dat.gz
	
# PROTEIN_ISOFORMS
	echo '##### Downloading PROTEIN_ISOFORMS'
	wget ftp://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/reference_proteomes/Eukaryota/UP000002311_559292_additional.fasta.gz
	
# PROTEIN_DOMAIN_SMART
	echo '##### Downloading PROTEIN_DOMAIN_SMART'
	wget http://smart.embl.de/smart/descriptions.pl -O all_protein_domains_smart.txt
	
# PROTEIN_DOMAIN_PFAM
	echo '##### Downloading PROTEIN_DOMAIN_PFAM'
	wget ftp://ftp.ebi.ac.uk/pub/databases/Pfam/current_release/Pfam-A.clans.tsv.gz
	
#
# [GENE_ONTONLOGY]
#
	
# GENE_ONTOLOGY_DEFINITION
	echo '##### Downloading GENE_ONTOLOGY_DEFINITION'
	wget http://purl.obolibrary.org/obo/go/go-basic.obo
	
# GENE_ONTOLOGY_ANNOTATION
	echo '##### Downloading GENE_ONTOLOGY_ANNOTATION'
	wget ftp://ftp.ebi.ac.uk/pub/databases/GO/goa/YEAST/gene_association.goa_yeast.gz
	
#
# [KEGG_PATHWAY]
#
	
# KEGG_PATHWAY_DEFINITION
	echo '##### Downloading KEGG_PATHWAY_DEFINITION'
	wget http://rest.kegg.jp/list/pathway/sce -O yeast_kegg_pathway_definitions.txt
	
# KEGG_PATHWAY_ANNOTATION
	echo '##### Downloading KEGG_PATHWAY_ANNOTATION'
	wget http://rest.kegg.jp/link/pathway/sce -O yeast_kegg_pathway_annotations.txt
	
#
# [REACTOME PATHWAY]
# 
	
# REACTOME_PATHWAY_DEFINITION
	echo '##### Downloading REACTOME_PATHWAY_DEFINITION'
	wget http://www.reactome.org/download/current/ReactomePathways.txt -O all_reactome_pathway_definitions.txt
	
# REACTOME_PATHWAY_ANNOTATION
	echo '##### Downloading REACTOME_PATHWAY_ANNOTATION'
	wget http://www.reactome.org/download/current/UniProt2Reactome.txt -O all_reactome_pathway_annotations.txt
	
#
# [INTERACTOME]
#
	
# See INTERACTOME section below in the "All three species" part 
	
	
#
# Unzip all compressed files
#
	gunzip $ROOT_DATA_FOLDER/$YEAST_DATA_FOLDER/*.gz
	
	
	
# ------------------------------------------------------------------------
# ALL THREE SPECIES
# ------------------------------------------------------------------------
	
	echo '#######################################'
	echo '### INTERACTOME FOR ALL THREE SPECIES #'
	echo '#######################################'
	
#
# [INTERACTOME]
#
# INTERACTOME_DEFINITION
# INTERACTOME_NETWORK_DEFINITION
	
# The production of the two above required files requires to execute MoonGo scripts.
#
# Running MoonGO requires to obtain the TREMBL protein sequences in the genbank flat file format (see http://www.ncbi.nlm.nih.gov/Sitemap/samplerecord.html)
# This flat file must be present in the MoonGO script folder. To obtain that file, the 'update_uniprot.sh' script is run. This script
# download the global flat files and extract the information for classical species using the script 'extract_species_from_flat.pl'. All the flat files for
# the studied species will be produced. 
	
	echo '##### Updating TREMBL protein sequence flat files'
	cd $ROOT_DATA_FOLDER/script/protein
	./update_uniprot.sh
	
# INTERACTOME_DEFINITION
# INTERACTOME_NETWORK_DEFINITION
	
# To run MoonGO, first export the $SPECIES environment variable to indicate the desired species (here human):
# export $SPECIES=human
# Then copy the flat files obtain in the previous step to the MoonGO script folder.
	
	echo '##### Cleaning MoonGO folder'
	cd $ROOT_MOONGO_FOLDER
	rm -f *.flat
	rm -f *.psi*
	rm -f *.gr
	rm -f *.pairmap
	rm -f *.blastmap
	rm -f *.missing
	rm -f error_hq.log
	rm -rf .psi_tmp
	
	echo '##### Copying flat files to MoonGO folder'
	cd $ROOT_DATA_FOLDER
	/bin/cp -f human.flat $ROOT_MOONGO_FOLDER
	/bin/cp -f fly.flat $ROOT_MOONGO_FOLDER
	/bin/cp -f yeast.flat $ROOT_MOONGO_FOLDER
	 
# Finally run the psciquic_wrapper.pl script:
#
# psicquic_wrapper.pl -FdBv --binary --net-file $SPECIES.binary.gr --flat-file $SPECIES.flat -p $SPECIES.psi --psi2id $SPECIES.psi2id -s $SPECIES -y 0.95
# psicquic_wrapper.pl -FdBv --net-file $SPECIES.cocomplex.gr --flat-file $SPECIES.flat -p $SPECIES.psi --psi2id $SPECIES.psi2id -s $SPECIES -y 0.95
#
	
	echo '##### Building interactomes and PPI Networks'
	cd $ROOT_MOONGO_FOLDER
	export PATH=$CDHIT_FOLDER:$ROOT_MOONGO_FOLDER:$ROOT_CLASS_ANNOTATION_FOLDER:$PATH
	
	echo '##### -- Executing MoonGO for HUMAN'
	export SPECIES=human
	./psicquic_wrapper.pl -FdBv --binary --net-file $SPECIES.binary.gr --flat-file $SPECIES.flat -p $SPECIES.psi --psi2id $SPECIES.psi2id -s $SPECIES -y 0.95
	echo '##### -- Copying result to HUMAN data folder'
	/bin/cp -f $SPECIES.pairmap $ROOT_DATA_FOLDER/$HUMAN_DATA_FOLDER
	/bin/cp -f $SPECIES.binary.gr $ROOT_DATA_FOLDER/$HUMAN_DATA_FOLDER
	/bin/cp -f $SPECIES.binary.nr0.95.gr $ROOT_DATA_FOLDER/$HUMAN_DATA_FOLDER
	/bin/cp -f ${SPECIES}_0.95.blastmap $ROOT_DATA_FOLDER/$HUMAN_DATA_FOLDER
	
	echo '##### -- Executing MoonGO for DROME'
	export SPECIES=fly
	./psicquic_wrapper.pl -FdBv --binary --net-file $SPECIES.binary.gr --flat-file $SPECIES.flat -p $SPECIES.psi --psi2id $SPECIES.psi2id -s $SPECIES -y 0.95
	echo '##### -- Copying result to DROME data folder'
	/bin/cp -f $SPECIES.pairmap $ROOT_DATA_FOLDER/$DROSOPHILA_DATA_FOLDER
	/bin/cp -f $SPECIES.binary.gr $ROOT_DATA_FOLDER/$DROSOPHILA_DATA_FOLDER
	/bin/cp -f $SPECIES.binary.nr0.95.gr $ROOT_DATA_FOLDER/$DROSOPHILA_DATA_FOLDER
	/bin/cp -f ${SPECIES}_0.95.blastmap $ROOT_DATA_FOLDER/$DROSOPHILA_DATA_FOLDER
	
	echo '##### -- Executing MoonGO for YEAST'
	export SPECIES=yeast
	./psicquic_wrapper.pl -FdBv --binary --net-file $SPECIES.binary.gr --flat-file $SPECIES.flat -p $SPECIES.psi --psi2id $SPECIES.psi2id -s $SPECIES -y 0.95
	echo '##### -- Copying result to YEAST data folder'
	/bin/cp -f $SPECIES.pairmap $ROOT_DATA_FOLDER/$YEAST_DATA_FOLDER
	/bin/cp -f $SPECIES.binary.gr $ROOT_DATA_FOLDER/$YEAST_DATA_FOLDER
	/bin/cp -f $SPECIES.binary.nr0.95.gr $ROOT_DATA_FOLDER/$YEAST_DATA_FOLDER
	/bin/cp -f ${SPECIES}_0.95.blastmap $ROOT_DATA_FOLDER/$YEAST_DATA_FOLDER


# The file '$SPECIES.pairmap' (INTERACTOME_DEFINITION) and the file '$SPECIES.<method>.gr' (INTERACTOME_NETWORK_DEFINITION) will be produced 
# (<method> should be 'binary' or 'cocomplex').
# Note that the non-redundant version of the '$SPECIES.<method>.gr' file will be produced with name '$SPECIES.<method>.nr0.95.gr' 
# and may be used as INTERACTOME_NETWORK_DEFINITION. The mapping between proteins that were identified as redundant is
# given in the $SPECIES_0.95.blasmap file. This file can be used to insert new crossreferences as INTERACTOME_NETWORK_REDUNDANCY_DEFINITION.


}
#############################################################################################
#############################################################################################
#### END OF FIRST PART FUNCTION: download_and_prepare_data()
#############################################################################################
#############################################################################################



#############################################################################################
#############################################################################################
#### START OF SECOND PART FUNCTION: classify_and_annotate()
#############################################################################################
#############################################################################################
classify_and_annotate(){


# INTERACTOME_NETWORK_PARTITION_DEFINITION

# The production of the partition file requires to execute OCG on a connex graph. So, the $SPECIES.<method>.gr or the $SPECIES.<method>.nr0.95.gr produced above
# must be opened in cytoscape to remove all non connex groups before running OCG (see tagc-ocg project on sourcesup). The resulting file
# will be named $SPECIES.<method>.connected.gr or the $SPECIES.<method>.nr0.95.connected.gr
# Note : <method> should be 'binary' or 'cocomplex' and the <species_folder> is $HUMAN_DATA_FOLDER, $DROSOPHILA_DATA_FOLDER or $YEAST_DATA_FOLDER

	echo '##### Partitioning PPI Networks'
	cd $ROOT_DATA_FOLDER

	echo '##### -- Executing OCG for HUMAN'
	export SPECIES=human
	/home/rainet/Software/ocg/OCG $ROOT_DATA_FOLDER/$HUMAN_DATA_FOLDER/$SPECIES.binary.nr0.95.connected.gr > $ROOT_DATA_FOLDER/$HUMAN_DATA_FOLDER/$SPECIES.binary.nr0.95.connected.clas
	
	echo '##### -- Executing OCG for DROME'
	export SPECIES=fly
	/home/rainet/Software/ocg/OCG $ROOT_DATA_FOLDER/$DROSOPHILA_DATA_FOLDER/$SPECIES.binary.nr0.95.connected.gr > $ROOT_DATA_FOLDER/$DROSOPHILA_DATA_FOLDER/$SPECIES.binary.nr0.95.connected.clas
	
	echo '##### -- Executing OCG for YEAST'
	export SPECIES=yeast
	/home/rainet/Software/ocg/OCG $ROOT_DATA_FOLDER/$YEAST_DATA_FOLDER/$SPECIES.binary.nr0.95.connected.gr > $ROOT_DATA_FOLDER/$YEAST_DATA_FOLDER/$SPECIES.binary.nr0.95.connected.clas

# INTERACTOME_NETWORK_PARTITION_ANNOTATION

# The production of the annotated partition file requires to execute an annotation script produced by Benoit Robisson (see tagc-ocg project on sourcesup).
# This script needs two files : the GO definition file (go-basic.obo downloaded in a previous step) and the annotation file (gene_association.goa... 
# downloaded in a previous step also). The go-basic.obo file must be copied to the 'annotation_files' folder located in the script folder.
# the gene_association.goa... file must be copied in the annotation_files/QuickGO folder
# The association files (gene_association.goa_$SPECIES) must be copied to the folder 'annotation_files/QuickGO'.

	echo '##### Copying GO definition and annotations to ClassAnnotation tool folders'
	/bin/cp -f $ROOT_DATA_FOLDER/$HUMAN_DATA_FOLDER/go-basic.obo $ROOT_CLASS_ANNOTATION_FOLDER/annotation_files
	/bin/cp -f $ROOT_DATA_FOLDER/$HUMAN_DATA_FOLDER/gene_association.goa_human $ROOT_CLASS_ANNOTATION_FOLDER/annotation_files/QuickGO
	/bin/cp -f $ROOT_DATA_FOLDER/$DROSOPHILA_DATA_FOLDER/gene_association.goa_fly $ROOT_CLASS_ANNOTATION_FOLDER/annotation_files/QuickGO
	/bin/cp -f $ROOT_DATA_FOLDER/$YEAST_DATA_FOLDER/gene_association.goa_yeast $ROOT_CLASS_ANNOTATION_FOLDER/annotation_files/QuickGO

# Once done, execute the annotation script (WARNING : DO NOT forget the -s option):

	echo '##### Annotating PPI Networks modules'
	cd $ROOT_CLASS_ANNOTATION_FOLDER 
	
	echo '##### -- Executing ClassAnnotation for HUMAN'
	export SPECIES=human
	python annotate.py -s $SPECIES $ROOT_DATA_FOLDER/$HUMAN_DATA_FOLDER/$SPECIES.binary.nr0.95.connected.clas > $ROOT_DATA_FOLDER/$HUMAN_DATA_FOLDER/$SPECIES.binary.nr0.95.connected.fm
	
	echo '##### -- Executing ClassAnnotation for DROME'
	export SPECIES=fly
	python annotate.py -s $SPECIES $ROOT_DATA_FOLDER/$DROSOPHILA_DATA_FOLDER/$SPECIES.binary.nr0.95.connected.clas > $ROOT_DATA_FOLDER/$DROSOPHILA_DATA_FOLDER/$SPECIES.binary.nr0.95.connected.fm
	
	echo '##### -- Executing ClassAnnotation for YEAST'
	export SPECIES=yeast
	python annotate.py -s $SPECIES $ROOT_DATA_FOLDER/$YEAST_DATA_FOLDER/$SPECIES.binary.nr0.95.connected.clas > $ROOT_DATA_FOLDER/$YEAST_DATA_FOLDER/$SPECIES.binary.nr0.95.connected.fm

}
#############################################################################################
#############################################################################################
#### END OF SECOND PART FUNCTION: classify_and_annotate()
#############################################################################################
#############################################################################################

 
# function to display menus
show_menus() {
	clear
	echo "~~~~~~~~~~~~~~~~~~~~~"	
	echo " M A I N - M E N U"
	echo "~~~~~~~~~~~~~~~~~~~~~"
	echo "1. Download and prepare data"
	echo "2. Classify and annotate connected PPI Networks"
	echo "3. Exit"
}
# read input from the keyboard and take a action
# invoke the one() when the user select 1 from the menu option.
# invoke the two() when the user select 2 from the menu option.
# Exit when user the user select 3 form the menu option.
read_options(){
	local choice
	read -p "Enter choice [ 1 - 3] " choice
	case $choice in
		1) download_and_prepare_data ;;
		2) classify_and_annotate ;;
		3) exit 0;;
		*) echo -e "${RED}Error...${STD}" && sleep 2
	esac
}
# Pause at the end of an option
pause(){
  read -p "######Process finished : press any key to continue...######" fackEnterKey
}
 
# ----------------------------------------------
# Trap CTRL+C, CTRL+Z and quit singles
# ----------------------------------------------
trap '' SIGINT SIGQUIT SIGTSTP
 
# -----------------------------------
# Main logic
# ------------------------------------
#while true
#do
	show_menus
	read_options
#done
