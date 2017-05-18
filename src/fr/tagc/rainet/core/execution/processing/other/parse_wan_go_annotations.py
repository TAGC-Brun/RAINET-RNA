
import re
import argparse

from fr.tagc.rainet.core.util.exception.RainetException import RainetException
from fr.tagc.rainet.core.util.log.Logger import Logger
from fr.tagc.rainet.core.util.time.Timer import Timer

from goatools.base import download_go_basic_obo
from goatools.obo_parser import GODag
from goatools.go_search import GoSearch


#===============================================================================
# Started 10-May-2017
# Diogo Ribeiro
DESC_COMMENT = "Script to process Wan 2015 complexes GO term annotation."
SCRIPT_NAME = "parse_wan_go_annotations.py"
#===============================================================================

#===============================================================================
# General plan:
# 1) Use goatools library to load Gene Ontology
# 2) Read Wan complex GO term annotations, and complexes we want to analyse
# 3) Produce frequencies of GO terms at a certain depth-level on our wanted complexes
#===============================================================================

#===============================================================================
# Processing notes:
# 1) For each GO term annotation in a complex, get its parents at the given depth-level
# 2)
#===============================================================================


# Global variables
theParents = set() # global variable storing parents of a term

# Constants
UNANNOTATED_CATEGORY = "'unannotated'"
PROPORTION_WAN_ANNOTATIONS = "proportion_wan_go_term_annotations"

# #
# Read GO term information (go basic file) using goatools
def read_go_basic(): 
    
    oboFile = download_go_basic_obo()
     
    obodag = GODag( oboFile)

    # Format example
    # ['name', 'level', 'is_obsolete', 'namespace', 'id', 'depth', 'parents', 'children', '_parents', 'alt_ids']
    # name:secondary active monocarboxylate transmembrane transporter activity
    # level:5
    # is_obsolete:False
    # namespace:molecular_function
    # id:GO:0015355
    # depth:9
    # parents: 2 items (more info)
    # children: 0 items
    # alt_ids: 0 items, GOTerm('GO:0042879'):

    return obodag
  
     
# #
# Read custom-made Wan 2015 GO term association file
def read_wan_go_term_annotation_file( inputFile, obodag):

    ##########################################
    # Read file with Wan complexes and their GO term annotations
    ##########################################
    # Example format
    # 507     plasma membrane (GO:0005886)
    # 296     isomerase activity (GO:0016853), cellular nitrogen compound metabolic process (GO:0034641)
    # Note: all types of ontologies (BP, CC, MF) were used
    # Note: one line per complex
    
    complexGOA = {} # key -> complexID, val -> list of associated GO terms
    
    # Read input file
    with open( inputFile, "r") as inFile:
        
        for line in inFile:
            line = line.strip()
            spl = line.split("\t")
    
            # Wan complex ID
            complexID = spl[0]
    
            # Each complex can have one or several
            # Note: sometimes terms themselves have a ",", therefore splitting by comma has to be done carefully        
            goTerms = spl[1].split(",")
                
            # Get all GO term IDs        
            goIDs = set()
            for term in goTerms:
                id = re.search( ".*(GO:[0-9]*).*", term)
                if id != None:
                    goIDs.add( id.group(1) )
            
            if complexID not in complexGOA:
                complexGOA[ complexID] = set()
            else:
                print "ERROR: duplicate complexID", complexID
                break
    
            complexGOA[ complexID] = goIDs

    print "Number of complexes read in annotation file: ", len( complexGOA)

    return complexGOA


# #
# recursively return all the parents of a certain go term, with a certain depth level
# Stores results into a previously defined global variable
def _get_parents( obodag, goTerm, wantedLevel):
    
    parents = obodag[ goTerm]._parents

    if obodag[ goTerm].level == wantedLevel:
        theParents.add( goTerm)
    else:
        for parent in parents:
            _get_parents(obodag, parent, wantedLevel)
    

# # 
# For each complex, loop its annotations and store only GO terms with a certain depth level (by checking their parents, if necessary)
def get_wanted_level_go_terms( complexGOA, obodag, wantedNamespace, wantedLevel):

    global theParents

    complexRoot = {} # Key -> complexID, val -> wanted level (parent) go terms 

    for complexID in complexGOA:
        for goTerm in complexGOA[ complexID]:
            
            assert len( goTerm) > 9
            
            goInfo = obodag[ goTerm]
     
            # Filter for wanted namespace
            if goInfo.namespace == wantedNamespace:
     
                if goInfo.level > wantedLevel:
                    # Search for parent of level 1
    
                    theParents = set() # important to reset global variable at each term
                    
                    # recursive function that modifies global variable to return all parents of a term
                    _get_parents( obodag, goTerm, wantedLevel)
                    
                    assert len( theParents) > 0
                        
                    # for each complex, store base level go terms             
                    complexRoot[ complexID] = theParents
                    
                else:
                    # If current GO term is of the wanted depth level, simply keep it
                    complexRoot[ complexID] = set([goTerm])

        # initialise complexes without annotation on the wantedNamespace
        if complexID not in complexRoot:
            complexRoot[ complexID] = set()

#     print "Number of complexes in annotation file: %s" % ( len( complexRoot))

    return complexRoot

# #
# Read list of complexes we want to analyse
def read_wanted_complex_list( inputListWantedComplexes):
    
    wantedComplexList = set()
    
    with open( inputListWantedComplexes, "r") as inFile:
        for line in inFile:
            wantedComplexList.add( line.strip())

    print "Number of wanted complexes: %s" % ( len( wantedComplexList))

    return wantedComplexList


# #
# Analyse frequency of terms per complex
def term_per_complex_analysis( complexRoot, wantedComplexList, obodag, outFile):
        
    termFreqs = {} # Key -> Term ID, val -> count
    termFreqs[ UNANNOTATED_CATEGORY] = 0

    for complex in complexRoot:
        
        # Filter out complexes which we do not want to analyse
        if complex not in wantedComplexList:
            continue
        
        terms = complexRoot[ complex]
        
        if len( terms) == 0:
            termFreqs[ UNANNOTATED_CATEGORY] += 1
        
        for term in terms:

            if term not in termFreqs:
                termFreqs[ term] = 0

            # TODO: check this
            termFreqs[ term] += (1.0 / len( terms) )

    # Note: the Wan annotation file only contains complexes with at least one annotation.    
    # Add unannotated category, containing also complexes in the wantedComplexList that are not in the Wan annotation file

    for complex in wantedComplexList:
        if complex not in complexRoot:
            termFreqs[ UNANNOTATED_CATEGORY] += 1

    summComplexes = sum( termFreqs.values())

    print "Number of summary terms: %s" % ( len( termFreqs))
    print "Number of complexes with terms (incl. unannotated): %s" % ( summComplexes)

    outFile.write("GO_ID\tGO_NAME\tPROPORTION\n")

    for term in termFreqs:
        if term != UNANNOTATED_CATEGORY:
            outFile.write( "%s\t%s\t%s\n" % ( term, obodag[ term].name, termFreqs[ term] / float( summComplexes) ) )
        else:
            outFile.write( "%s\t%s\t%s\n" % ( term, UNANNOTATED_CATEGORY, termFreqs[ term] / float( summComplexes) ) )
            

    outFile.close()


# # #
# # Returns list of terms with the wanted depth level and namespace
# def read_terms_with_depth_level(obodag, wantedLevel, wantedNamespace):
# 
#     termsWithLevel = set() # set of GO terms on the wanted depth level
#     for term in obodag:
#         
#         if obodag[term].level == wantedLevel:
#             if obodag[term].namespace == wantedNamespace:
#                 termsWithLevel.add( term)
#                 print term, obodag[term].alt_ids
# 
#     print termsWithLevel
#     # the list includes GO IDs that are alternative IDs, therefore, it can be larger than the actual 'usual' list of level 1 terms (i.e. Quick GO)
#     
#     return termsWithLevel

#####################
####################
# Client
####################
#####################

if __name__ == "__main__":

    try:
    
        # Start chrono
        Timer.get_instance().start_chrono()
        print "STARTING " + SCRIPT_NAME
        
        #===============================================================================
        # Get input arguments
        #===============================================================================
        parser = argparse.ArgumentParser(description= DESC_COMMENT) 
    
        # positional args
        parser.add_argument('wanAnnotations', metavar='wanAnnotations', type=str,
                             help='Wan 2015 complex GO term annotation file (custom-made TSV). E.g. 818     lipid particle (GO:0005811), cytoskeleton organization (GO:0007010)')
        parser.add_argument('wantedComplexes', metavar='wantedComplexes', type=str,
                             help='List of Wan complexes which we want to analyse')
        parser.add_argument('--wantedNamespace', metavar='wantedNamespace', type=str, default = "biological_process",
                             help='Ontology type we want to analyse (biological_process, cellular_components, molecular_function.')
        parser.add_argument('--wantedLevel', metavar='wantedLevel', type=int, default = 1,
                             help='Gene Ontology depth-level we want to analyse (0,1,2,3,4 etc).')
           
        #gets the arguments
        args = parser.parse_args( ) 
    
        #===============================================================================
        # Run analysis / processing
        #===============================================================================

        outFile = open( PROPORTION_WAN_ANNOTATIONS + "_" + args.wantedNamespace + "_" + str( args.wantedLevel) + ".txt", "w")

        obodag = read_go_basic()
#         read_terms_with_depth_level( obodag, args.wantedLevel, args.wantedNamespace)
        complexGOA = read_wan_go_term_annotation_file( args.wanAnnotations, obodag)
        complexRoot = get_wanted_level_go_terms( complexGOA, obodag, args.wantedNamespace, args.wantedLevel)
        wantedComplexList = read_wanted_complex_list( args.wantedComplexes)
        term_per_complex_analysis( complexRoot, wantedComplexList, obodag, outFile)


        # Stop the chrono      
        Timer.get_instance().stop_chrono( "FINISHED " + SCRIPT_NAME )

    # Use RainetException to catch errors
    except RainetException as rainet:
        Logger.get_instance().error( "Error during execution of %s. Aborting :\n" % SCRIPT_NAME + rainet.to_string())

