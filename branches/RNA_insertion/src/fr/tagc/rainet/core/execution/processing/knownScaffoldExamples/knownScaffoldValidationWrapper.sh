
script='python2.7 /home/diogo/workspace/tagc-rainet-RNA/src/fr/tagc/rainet/core/execution/processing/knownScaffoldExamples/KnownScaffoldValidation.py'

catrapid_folder='/home/diogo/Documents/RAINET_data/catRAPID/webserver_results/'

#npinter_file='/home/diogo/Documents/RAINET_data/NPInter/golden_set_NPInter[v3.0].txt'
npinter_file='/home/diogo/Documents/RAINET_data/NPInter/interaction_NPInter[v3.0].txt'

rna_list_folder='/home/diogo/Documents/RAINET_data/TAGC/rainetDatabase/results/knownScaffoldValidation/wantedRNAFiles/'

rainet_db='/home/diogo/Documents/RAINET_data/TAGC/rainetDatabase/db_backup/RNA/rainet2016-03-17.human_expression.sqlite'

out_folder='/home/diogo/Documents/RAINET_data/TAGC/rainetDatabase/results/knownScaffoldValidation/scoreColumn'

topProportion='1'

scoreColumn=1


#### XIST MOUSE ####

$script $catrapid_folder/all_interactions.85891.XIST_MOUSE_RBP_DBP_DISORDER.txt /home/diogo/Documents/RAINET_data/TAGC/rainetDatabase/results/knownScaffoldValidation/manualLists/XIST_paper_table_s1_tab2_uniprot.txt fakeFile $rainet_db $out_folder$scoreColumn --npinter 0 --columnForPlot $scoreColumn
# $script $catrapid_folder/all_interactions.85891.XIST_MOUSE_RBP_DBP_DISORDER.txt /home/diogo/Documents/RAINET_data/TAGC/rainetDatabase/results/knownScaffoldValidation/manualLists/XIST_Tartaglia_paper_table_s1c.txt fakeFile $rainet_db $out_folder$scoreColumn --npinter 0 --columnForPlot $scoreColumn

#### SAMMSON ####

list_file="SAMMSON_keywords.tsv"

$script $catrapid_folder/ENST00000483525_SAMMSON.tsv /home/diogo/Documents/RAINET_data/TAGC/rainetDatabase/results/knownScaffoldValidation/manualLists/SAMMSON_paper_table_s2_tab2_uniprot.txt fakeFile $rainet_db $out_folder$scoreColumn --npinter 0 --columnForPlot $scoreColumn --searchSpace /home/diogo/Documents/RAINET_data/TAGC/rainetDatabase/results/knownScaffoldValidation/manualLists/SAMMSON_paper_table_s2_all_ids.txt

#exit

#### NORAD ####

list_file="NORAD_keywords.tsv"

$script $catrapid_folder/ENST00000565493_NORAD.tsv /home/diogo/Documents/RAINET_data/TAGC/rainetDatabase/results/knownScaffoldValidation/manualLists/NORAD_paper_table_s1.txt fakeFile $rainet_db $out_folder$scoreColumn --npinter 0 --columnForPlot $scoreColumn

#### NEAT1 ####

list_file="NEAT1_keywords.tsv"

$script $catrapid_folder/ENST00000501122_NEAT1_001.tsv $npinter_file $rna_list_folder$list_file $rainet_db $out_folder$scoreColumn  --columnForPlot $scoreColumn	
$script $catrapid_folder/ENST00000499732_NEAT1_002.tsv $npinter_file $rna_list_folder$list_file $rainet_db $out_folder$scoreColumn  --columnForPlot $scoreColumn

#### HOTAIR ####

list_file="HOTAIR_keywords.tsv"

$script $catrapid_folder/ENST00000424518_HOTAIR_001.tsv $npinter_file $rna_list_folder$list_file $rainet_db $out_folder$scoreColumn  --columnForPlot $scoreColumn
$script $catrapid_folder/ENST00000455246_HOTAIR_002.tsv $npinter_file $rna_list_folder$list_file $rainet_db $out_folder$scoreColumn  --columnForPlot $scoreColumn

#### MALAT1 ####

list_file="MALAT1_keywords.tsv"

$script $catrapid_folder/ENST00000534336_MALAT1_001.tsv $npinter_file $rna_list_folder$list_file $rainet_db $out_folder$scoreColumn  --columnForPlot $scoreColumn

#### MEG3 ####

list_file="MEG3_keywords.tsv"

$script $catrapid_folder/ENST00000423456_MEG3_002.tsv $npinter_file $rna_list_folder$list_file $rainet_db $out_folder$scoreColumn  --columnForPlot $scoreColumn
$script $catrapid_folder/ENST00000522771_MEG3_015.tsv $npinter_file $rna_list_folder$list_file $rainet_db $out_folder$scoreColumn  --columnForPlot $scoreColumn

#### KCNQ1OT1 ####

list_file="KCNQ1OT1_keywords.tsv"

$script $catrapid_folder/ENST00000597346_KCNQ1OT1_001.tsv $npinter_file $rna_list_folder$list_file $rainet_db $out_folder$scoreColumn  --columnForPlot $scoreColumn

#### XIST ####

list_file="XIST_keywords.tsv"

$script $catrapid_folder/ENST00000429829_XIST_001.tsv $npinter_file $rna_list_folder$list_file $rainet_db $out_folder$scoreColumn  --columnForPlot $scoreColumn

