#!/usr/bin/env python
##############################################################################
# Copyright 2016  Olexiouk Volodimir,Menschaert Gerben                       #
#                                                                            #
# Licensed under the Apache License, Version 2.0 (the "License");            #
# you may not use this file except in compliance with the License.           #
# You may obtain a copy of the License at                                    #
#                                                                            #
#  http://www.apache.org/licenses/LICENSE-2.0                                #
#                                                                            #
# Unless required by applicable law or agreed to in writing, software        #
# distributed under the License is distributed on an "AS IS" BASIS,          #
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.   #
# See the License for the specific language governing permissions and        #
# limitations under the License.                                             #
#                                                                            #
##############################################################################

from __future__ import division

__author__ = 'Volodimir Olexiouk'

import time
import pysam
import pysam.ctabixproxies
import proBAM_input
import proBAM_ENSEMBL
import sys
import os
import proBAM_proBED
import proBAM_IDparser
import argparse
from proBAM_coref import *
try:
    from collections import OrderedDict
except ImportError:
    from ordereddict import OrderedDict

######################################
###         DEPENDENCIES           ###
######################################
'''
cogent
pysam
mySQLdb
SQLAlchemy
lxml
numpy
matplotlib
pyteomics
'''
#
# Parse command line arguments
#
def get_parser():
    '''
    :return: parser
    '''
    parser = argparse.ArgumentParser(description=('proBAMconvert version 1.0.0'))
    parser.add_argument('--name','-N',help='name of the project (will be determine how the output file is called',
                        required=True)
    parser.add_argument('--mismatches', '-M', help='numpber of mismatches allowed during mapping',
                        required=False,default=0,choices=['0','1','2','3','4','5'],dest='allowed_mismatches')
    parser.add_argument('--version', '-V', help='ENSEMBL version to be used',
                        required=False,default=87,choices=['74','75','76','77','78','79','80','81','82','83','84','85','86','87'],
                        dest='database_v')
    parser.add_argument('--database', '-D', help='Which database has to be used (currently only ENSEMBL is available',
                        default="ENSEMBL",choices=['ENSEMBL'])
    parser.add_argument('--species', '-S', help='species to be used',
                        required=True,choices=['homo_sapiens','mus_musculus','drosophila_melanogaster','danio_rerio'])
    parser.add_argument('--file', '-F', help='location of the psm file to be processed',
                        required=True,dest='psm_file')
    parser.add_argument('--directory', '-C', help='location where to output files should be stored',
                        default=os.getcwd())
    parser.add_argument('--rm_duplicates', '-R', help='Whether duplicates should be removes',
                        default='N',choices=['Y','N'])
    parser.add_argument('--three_frame_translation', '-T', help='translate transcript sequences in 3 frames',
                        default='N',choices=['Y','N'])
    parser.add_argument('--decoy_annotation', '-E', help='annotation for DECOY PSM',
                        default='REV_,DECOY_,_REVERSED,REVERSED_,_DECOY')
    parser.add_argument('--sorting_order', '-O', help='sorting order of the SAM file',
                        default='unknown',choices=['unknown','unsorted','queryname','coordinate'])
    parser.add_argument('--pre_picked_annotation', '-P', help='Which/How annotation should be identified',
                        default='First',choices=['First','Ensembl_tr', 'Ensembl_pr','UniProt_ACC','UniProt_Entry',
                    'RefSeq','all'])
    parser.add_argument('--include_unmapped', '-U', help='Whether unmapped psm should be included in the output',
                        default='Y',choices=['Y','N'])
    parser.add_argument('--conversion_mode', '-X', help='which ouput format should be generated',
                        default='proBAM_psm',choices=['proBAM_psm','proBAM_peptide','proBAM_peptide_mod','proBED'])
    parser.add_argument('--comments', '-Y', help='add a comment to the file',
                        default='')
    return parser

#
# Command line variable input
#
def get_input_variables():

    ###############################
    # GLOBAL VARIABLE DECLARATION #
    ###############################

    global name
    global allowed_mismatches
    global database_v
    global database
    global species
    global psm_file
    global directory
    global decoy_annotation
    global sorting_order
    global rm_duplicates
    global three_frame_translation
    global command_line
    global comments
    global conversion_mode
    global pre_picked_annotation
    global include_unmapped
    global version

    ######################################
    ### VARIABLE AND INPUT DECLARATION ###
    ######################################
    version='1.0.0'
    decoy_annotation=decoy_annotation.split(',')

    command_line = "python proBAM.py --name " + str(name) + " --mismatches " + str(
                   allowed_mismatches) + " --version " + str(database_v) \
                   + " --database " + str(database) + " --species " + str(species) + " --file " + str(psm_file) + \
                   " --directory " + str(directory) + " --rm_duplicates " + str(rm_duplicates) + \
                   " --tri_frame_translation " + str(three_frame_translation)+\
                   " --conversion_mode "+str(conversion_mode)

    # ouput variables
    print(  "psm file:                                      " + str(psm_file) +"\n"+
            "directory:                                     " + str(directory) +"\n"+
            "database                                       " + str(database) +"\n"+
            "database version:                              " + str(database_v) +"\n"+
            "species:                                       " + str(species) +"\n"+
            "allowed mismatches:                            " + str(allowed_mismatches)+"\n"+
            "three_frame_translation:                       " + str(three_frame_translation)+"\n"+
            "decoy annotations:                             " + str(decoy_annotation)+"\n"+
            "sorting order:                                 " + str(sorting_order)+"\n"+
            "pre picked annotation                          " + str(pre_picked_annotation)+"\n"+
            "conversion_mode:                               " + str(conversion_mode))


'''
###############################
# NON COMMAND LINE ARGUMENTS  #
# FOR TESTING PURPOSES        #
###############################

directory="/home/vladie/Desktop/proBAMconvert/output/"
psm_file="/home/vladie/Desktop/proBAMconvert/PXD000652.mzid"
species="homo_sapiens"
database='ENSEMBL'
database_v=87
# TODO Let users specify used the decoy annotation
decoy_annotation=['REV_','DECOY_','_REVERSED','REVERSED_','_DECOY']
allowed_mismatches=0
version='1.0'
# can be unknown,unsorted, queryname or coordinate, can be specified by user
sorting_order='unknown'
name='test4'
three_frame_translation='N'
rm_duplicates="Y"
probed='N'
comments=''
include_unmapped='N'
pre_picked_annotation="all"
conversion_mode='proBAM_psm'

command_line= "python proBAM.py --name "+str(name)+" --mismatches "+str(allowed_mismatches)+" --version "+str(database_v)\
              +" --database "+str(database)+" --species "+str(species)+" --file "+str(psm_file)+\
              " --directory "+str(directory)+" --rm_duplicates "+str(rm_duplicates)+\
              " --tri_frame_translation "+str(three_frame_translation+
              "--pre_picked_annotation "+str(pre_picked_annotation)+" --conversion_mode "+str(conversion_mode))

# ouput variables
print(  "psm file:                                      " + str(psm_file) +"\n"+
        "directory:                                     " + str(directory) +"\n"+
        "database                                       " + str(database) +"\n"+
        "database version:                              " + str(database_v) +"\n"+
        "species:                                       " + str(species) +"\n"+
        "allowed mismatches:                            " + str(allowed_mismatches)+"\n"+
        "three_frame_translation:                       " + str(three_frame_translation)+"\n"+
        "remove duplicates:                             " + str(rm_duplicates)+"\n"+
        "pre picked annotation:                         " + str(pre_picked_annotation)+"\n"+
        "include_unmapped:                              " + str(include_unmapped))
'''
#######################
### GETTERS/SETTERS ###
#######################

#
# Function to write temp_result to file
#
def write_psm(temp_result,file):
    '''
    :param temp_result: sam results
    :param file: sam file
    :return: write to sam file IO
    '''
    for tab in temp_result:
        file.write(str(tab)+'\t')
    file.write('\n')

#
# Open file to write to
#
def open_sam_file(directory,name):
    '''
    :param directory: working directory
    :param name: output file name
    :return: return output file IO
    '''
    file=open(directory+name+'.sam',"w")
    return file


######################
### MAIN FUNCTIONS ###
######################

#
# Convert PSM to SAM
#
def PSM2SAM(psm_hash,transcript_hash,exon_hash,decoy_annotation,allowed_mismatches,file,rm_duplicates,
            three_frame_translation,psm_file,id_map,gui):
    '''
    :param psm_hash: dictionairy of psm files
    :param transcript_hash: dictionairy of transcripts
    :param exon_hash: dictionairy of exons
    :param decoy_annotation: decoy annotation list
    :param allowed_mismatches: number of allowed mismatches
    :param file: sam file
    :return: sam file IO
    '''
    print("Commencing generation of SAM file")
    # psm_hash.reset()
    if rm_duplicates == "Y":
        dup = {}
    enzyme=get_enzyme(psm_file)
    enzyme_specificity=get_enzyme_specificity(psm_file)
    total_psms=len(psm_hash)
    current_psm=1
    percentage=5
    to_write=[]
    print('total psm_2B_processed:', len(psm_hash))
    for psm in psm_hash:
        # track progress
        current_psm+=1
        if float(current_psm*100)/float(total_psms) > percentage:
            print("%s%%" % percentage)
            percentage += 5
        # update window if in GUI
        if gui != None:
            gui.update()
        # convert unmapped PSMs with their own converter
        if 'search_hit' not in psm:
            continue
        else:
            for row in psm['search_hit']:
                for p in range(0, len(row['proteins'])):
                    decoy=0
                    # convert decoys with decoy-specific convertor
                    if 'DECOY_' in row['proteins'][p]['protein'].upper():
                        key = row['proteins'][p]['protein']
                        decoy=1
                        temp_result= decoy_PSM_to_SAM(psm,row,key,enzyme,enzyme_specificity)
                        if rm_duplicates=="Y":
                            dup_key= str(temp_result[0])+"_"+str(temp_result[9])+"_"+str(temp_result[2])\
                                     +"_"+str(temp_result[3])
                            if dup_key not in dup:
                                dup[dup_key]=1
                                to_write.append(temp_result)
                        else:
                            to_write.append(temp_result)

                    if decoy == 0:
                        key=row['proteins'][p]['protein']
                        # Filter out PSM where transcript sequences were not found/ non-existent
                        if (key not in id_map):
                            write_psm(unannotated_PSM_to_SAM(psm,row,decoy,key,enzyme,enzyme_specificity),file)
                        # transcript not on an canonical transcript
                        # TODO do this nicer by fetching canonical chr
                        else:
                            is_hit=0
                            for transcript_id in id_map[key]:
                                if (transcript_id not in transcript_hash or
                                            len(transcript_hash[transcript_id]['chr']) > 4):
                                    continue
                                elif three_frame_translation=="Y":
                                    temp_hit=map_peptide_to_protein_3frame(row['peptide'],
                                                                              transcript_hash[transcript_id]['transcript_seq'],
                                                                              allowed_mismatches,
                                                                              transcript_hash[transcript_id]['strand'])
                                    protein_hit=temp_hit[0]
                                    pre_post_aa=temp_hit[1]
                                else:
                                    temp_hit=map_peptide_to_protein(row['peptide'],transcript_hash[transcript_id]['protein_seq']
                                                                       ,allowed_mismatches)
                                    protein_hit=temp_hit[0]
                                    pre_post_aa=temp_hit[1]
                                if len(protein_hit)==0:
                                    continue
                                else:
                                    # map peptide on protein and retrieve hit position, iterate over all hits
                                    for phit in protein_hit:
                                        start_time = time.time()
                                        is_hit=1
                                        temp_result=[None]*33
                                        #
                                        # Mandatory columns adapted from SAM/BAM format
                                        #
                                        #QNAME
                                        temp_result[0]=psm['spectrum']
                                        #FLAG
                                        temp_result[1]=calculate_FLAG(transcript_hash[transcript_id]['strand'],row['hit_rank'],
                                                                      decoy)
                                        #RNAME
                                        temp_result[2]='chr'+str(transcript_hash[transcript_id]['chr'])
                                        #POS
                                        pos_and_exons=calculate_genome_position(phit[0],
                                                                 transcript_hash[transcript_id]['strand'],
                                                                 transcript_hash[transcript_id]['5UTR_offset'],
                                                                 transcript_hash[transcript_id]['start_exon_rank'],
                                                                 row['peptide'],
                                                                 exon_hash[transcript_hash[transcript_id]['transcript_id']],
                                                                 transcript_hash[transcript_id]['chr'],
                                                                 three_frame_translation,
                                                                transcript_hash[transcript_id]['shift'])
                                        temp_result[3]=pos_and_exons[0]
                                        #MAPQ
                                        temp_result[4]=255
                                        #CIGAR
                                        cigar = compute_cigar(temp_result[3], pos_and_exons[1],
                                                              transcript_hash[transcript_id]['strand'], row['peptide'])
                                        temp_result[5]=cigar[0]
                                        #RNEXT
                                        temp_result[6]='*'
                                        #PNEXT
                                        temp_result[7]=0
                                        #TLEN
                                        temp_result[8]=0
                                        #SEQ
                                        if three_frame_translation=='Y':
                                            phit_loc=phit[0]
                                        else:
                                            phit_loc=phit[0]*3
                                        if int(transcript_hash[transcript_id]['strand'])==1:
                                            temp_result[9]=str(transcript_hash[transcript_id]['transcript_seq']\
                                                       [phit_loc:(phit_loc+(len(row['peptide'])*3))])
                                        else:
                                            temp_result[9]=reverse_complement(str(transcript_hash[transcript_id]['transcript_seq']\
                                                       [phit_loc:(phit_loc+(len(row['peptide'])*3))]))
                                        #QUAL
                                        temp_result[10]='*'
                                        #
                                        #Mandatory proteomics specific columns added to the proBam format
                                        #
                                        #NH: number of genomic location the peptide mapping to
                                        temp_result[11]='NH:i:*'
                                        #XA: Whether the peptide is well annotated
                                        temp_result[12]=create_XA(phit[1])
                                        #XB: Mass error (experimental - calculated)
                                        temp_result[13]="XB:f:"+str(row['massdiff'])
                                        #XC: Peptide charge
                                        temp_result[14]='XC:i:'+str(psm['assumed_charge'])
                                        #XE: enzyme used
                                        temp_result[15]="XE:i:"+str(enzyme)
                                        #XF: reading frame of the peptide
                                        temp_result[16]='XF:Z:'+cigar[1]
                                        #XG: Petide type
                                        temp_result[17]=create_XG(phit[1])
                                        #XI: peptide intensity
                                        temp_result[18]="XI:f:*"
                                        #XL: number of peptides the spectrum mapping to
                                        temp_result[19]='XL:i:*'
                                        #XM: Modification
                                        temp_result[20]='XM:Z:'+create_XM(row['modifications'])
                                        #XN: number of mis-cleavages
                                        if 'num_missed_cleaveges' in row:
                                            temp_result[21]='XN:i:'+str(row['num_missed_cleavages'])
                                        else:
                                            temp_result[21]='XN:i:0'
                                        #XO: uniqueness of peptide mapping
                                        #todo figure this one out
                                        temp_result[22]='XO:Z:*'
                                        #XP; peptide sequence
                                        temp_result[23]='XP:Z:'+row['peptide']
                                        #XQ: PSM-Qvalue
                                        temp_result[24]='XQ:f:'+str(row['search_score']['evalue'])
                                        #XR: reference peptide sequence
                                        temp_result[25]='XR:Z:'+translate_seq(temp_result[9],
                                                                              transcript_hash[transcript_id]['strand'])
                                        #XS: PSM score
                                        temp_result[26]="XS:f:"+str(row['search_score']['score'])
                                        #XT: non/semi/tryptic
                                        temp_result[27]="XT:i:"+str(enzyme_specificity)
                                        #XU: petide URL
                                        temp_result[28]="XU:Z:*"
                                        #YA: following 2AA:
                                        temp_result[29]="YA:Z:"+str(pre_post_aa[1])
                                        #YB: preceding 2AA
                                        temp_result[30]="YB:Z:"+str(pre_post_aa[0])
                                        #YP: protein accession ID from the original search
                                        temp_result[31]='YP:Z:'+str(key)
                                        # ZA additional field specifiying the transcript/protein id used for mapping
                                        temp_result[32] = "ZA:Z:" + str(transcript_id)
                                        # remove duplicates if rm_duplicates=Y
                                        if rm_duplicates=="Y":
                                            dup_key= str(temp_result[9])+"_"+\
                                                              str(str(temp_result[0])+"_"+temp_result[2])+"_"+str(temp_result[3])
                                            if dup_key not in dup:
                                                dup[dup_key]=1
                                                to_write.append(temp_result)
                                        else:
                                            to_write.append(temp_result)
                            if is_hit==0:
                                to_write.append(unannotated_PSM_to_SAM(psm, row, decoy, key, enzyme, enzyme_specificity))
    print(" ")
    print("Writing SAM-file")
    for line in to_write:
        write_psm(line,file)
    file.close()

#
# Function to convert unannotated PSMs to SAM
#
def unannotated_PSM_to_SAM(psm,row,decoy,key,enzyme,enzyme_specificity):
    '''
    :param psm: psm dictionairy
    :param row: unnanotated PSM row
    :param decoy: decoy boolean
    :param file: output file
    :return: sam of unnanotated PSM
    '''
    decoy=int(decoy)
    temp_result=[None]*33
    #
    # Mandatory columns adapted from SAM/BAM format
    #
    #QNAME
    temp_result[0]=psm['spectrum']
    #FLAG
    temp_result[1]='4'
    #RNAME
    temp_result[2]='*'
    #POS
    temp_result[3]=0
    #MAPQ
    temp_result[4]=255
    #CIGAR
    temp_result[5]='*'
    #RNEXT
    temp_result[6]='*'
    #PNEXT
    temp_result[7]=0
    #TLEN
    temp_result[8]=0
    #SEQ
    temp_result[9]='*'
    #QUAL
    temp_result[10]='*'
    #
    #Mandatory proteomics specific columns added to the proBam format
    #
    #NH: number of genomic location the peptide mapping to
    temp_result[11]='NH:i:*'
    #XA: Whether the peptide is well annotated
    temp_result[12]='XA:i:2'
    #XB: Mass error
    temp_result[13]="XB:f:"+str(row['massdiff'])
    #XC: Peptide charge
    temp_result[14]='XC:i:'+str(psm['assumed_charge'])
    #XE: enzyme
    temp_result[15]="XE:i:"+str(enzyme)
    #XF: Reading frame of the peptide
    temp_result[16]="XF:Z:*"
    #XG: Petide type
    if decoy==1:
        temp_result[17]="XG:Z:D"
    else:
        temp_result[17]="XG:Z:U"
    #XI: Peptide intensity
    temp_result[18]="XI:f:*"
    #XL: number of peptides the spectrum mapping to
    temp_result[19]='XL:i:*'
    #XM: Modification
    temp_result[20]='XM:Z:'+create_XM(row['modifications'])
    #XN: number of mis-cleavages
    if 'num_missed_cleavages' in row:
        temp_result[21]='XN:i:'+str(row['num_missed_cleavages'])
    else:
        temp_result[21]='XN:i:*'
    #XO: uniqueness
    temp_result[22]='XO:Z:*'
    #XP; peptide sequence
    temp_result[23]='XP:Z:'+row['peptide']
    # XQ: PSM-Qvalue
    temp_result[24] = 'XQ:f:' + str(row['search_score']['evalue'])
    #XR: reference peptide sequence
    temp_result[25]='XR:Z:*'
    # XS: PSM score
    temp_result[26] = "XS:f:" + str(row['search_score']['score'])
    #XT: non/semi/tryptic
    temp_result[27]="XT:i:"+str(enzyme_specificity)
    #XU
    temp_result[28]="XU:Z:*"
    #YA: 2 AA after
    temp_result[29]='YA:Z:*'
    #YB: 2 AA before
    temp_result[30]='YB:Z:*'
    #YP: protein accession id
    temp_result[31]="YP:Z:"+str(key)
    #ZA additional field specifiying the transcript/protein id used for mapping
    temp_result[32] = "ZA:Z:*"
    return temp_result

#
# Function to convert decoy PSM to SAM format
#

def decoy_PSM_to_SAM(psm,row,key,enzyme,enzyme_specificity):
    '''
    :param psm: psm dictionairy
    :param row: row where decoy found
    :param key: psm key
    :param transcript_hash: transcript dictionairy
    :param exon_hash: exon dictionairy
    :param allowed_mismatches: number of allowed mismatches
    :param file: output file
    :return: SAM of decoy PSM
    '''

    # LEGACY: map decoy to genome if map_decoy=Y
    return unannotated_PSM_to_SAM(psm, row, 1, key, enzyme, enzyme_specificity)
    '''
    temp_result=[None]*23

    if map_decoy=="Y":
        protein_hit=map_peptide_to_protein(row['peptide'][::-1],transcript_hash[id_map[key]]['protein_seq'],allowed_mismatches)[0]
        pre_post_aa=map_peptide_to_protein(row['peptide'][::-1],transcript_hash[id_map[key]]['protein_seq'],allowed_mismatches)[1]
    else:

    protein_hit=[]
    if len(protein_hit)==0:
         return unannotated_PSM_to_SAM(psm,row,1,key,enzyme,enzyme_specificity)
    else:
        # map peptide on protein and retrieve hit position, iterate over all hits
        for phit in protein_hit:
            temp_result=[None]*32
            #
            # Mandatory columns adapted from SAM/BAM format
            #
            #QNAME
            temp_result[0]=psm['spectrum']
            #FLAG
            temp_result[1]=calculate_FLAG(transcript_hash[id_map[key]]['strand'],row['hit_rank'],
                                          1)
            #RNAME
            temp_result[2]='chr'+str(transcript_hash[id_map[key]]['chr'])
            #POS
            temp_result[3]=calculate_genome_position(phit[0],
                                                     transcript_hash[id_map[key]]['strand'],
                                                     transcript_hash[id_map[key]]['5UTR_offset'],
                                                     transcript_hash[id_map[key]]['start_exon_rank'],
                                                     row['peptide'][::-1],
                                                     exon_hash[transcript_hash[id_map[key]]['transcript_id']],
                                                     transcript_hash[id_map[key]]['chr'],
                                                     three_frame_translation)
            #MAPQ
            temp_result[4]=255
            #CIGAR
            temp_result[5]=compute_cigar(temp_result[3],
                                         exon_hash[transcript_hash[id_map[key]]['transcript_id']],
                                         transcript_hash[id_map[key]]['strand'],row['peptide'])
            #RNEXT
            temp_result[6]='*'
            #PNEXT
            temp_result[7]=0
            #TLEN
            temp_result[8]=0
            #SEQ
            if int(transcript_hash[id_map[key]]['strand'])==1:
                temp_result[9]=str(transcript_hash[id_map[key]]['transcript_seq']\
                           [(phit[0]*3):((phit[0]*3)+(len(row['peptide'])*3))])
            else:
                temp_result[9]=reverse_complement(str(transcript_hash[id_map[key]]['transcript_seq']\
                           [(phit[0]*3):((phit[0]*3)+(len(row['peptide'])*3))]))
            #QUAL
            temp_result[10]='*'
            #
            #Mandatory proteomics specific columns added to the proBam format
            #
            #NH: number of genomic location the peptide mapping to
            temp_result[11]='NH:i:'+str(len(row['proteins'])+len(phit)-1)
            # todo figure this one out
            temp_result[12] = 'XO:z:*'
            # XL: number of peptides the spectrum mapping to
            temp_result[13] = 'XL:i:*'
            # XP; peptide sequence
            temp_result[14] = 'XP:Z:' + row['modified_peptide']
            # YP: protein accession ID from the original search
            temp_result[15] = 'YP:Z:' + str(key)
            # XF: reading frame of the peptide
            temp_result[16] = 'XF:Z:' + compute_cigar(temp_result[3],
                                                      exon_hash[transcript_hash[id_map[key]]['transcript_id']],
                                                      transcript_hash[id_map[key]]['strand'], row['peptide'])[1]
            # XI: peptide intensity
            temp_result[17] = "XI:f:*"
            # XB: Mass error (experimental - calculated)
            temp_result[18] = "XB:f:" + str(row['massdiff'])
            # XR: reference peptide sequence
            temp_result[19] = 'XR:Z:' + row['peptide']
            # YB: preceding 2AA
            temp_result[20] = "YB:Z:*"
            # YA: following 2AA:
            temp_result[21] = "YA:Z:*"
            # XS: PSM score
            temp_result[22] = "XS:f:" + str(row['search_score']['score'])
            # XQ: PSM-Qvalue
            temp_result[23] = 'XQ:f:' + str(row['search_score']['evalue'])
            #XC: Peptide charge
            temp_result[24]='XC:i:'+str(psm['assumed_charge'])
            #XA: Whether the peptide is well annotated
            temp_result[25]=create_XA(phit[1])
            #XM: Modification
            temp_result[26]='XM:Z:'+create_XM(row['modifications'])
            #XN: number of mis-cleavages
            if 'num_missed_cleaveges' in row.keys():
                temp_result[27]='XN:i:'+str(row['num_missed_cleavages'])
            else:
                temp_result[27]='XN:i:*'
            #XT: non/semi/tryptic
            temp_result[28]="XT:i:"+str(enzyme_specificity)
            #XE enzyme
            temp_result[29]="XE:i"+str(enzyme)
            #XG: Petide type
            temp_result[30]='XG:Z:D'
            #XU= url
            temp_result[31]="XU:Z:*"
            return temp_result
    '''
#
# Create SAM header
#
def create_SAM_header(file,version,database,sorting_order,database_v,species,command_line,psm_file,comments):
    '''
    :param file: output file
    :param version: proBAMconvert version
    :param database: database name
    :param sorting_order: SAM sorting order
    :param database_v: database version
    :return:
    '''
    print('Creating SAM header')
    header=[]
    header.append('@HD\tVN:'+str(version)+' SO:'+sorting_order)
    if database.upper()=="ENSEMBL":
        SQ=proBAM_ENSEMBL.create_SQ_header(database_v,species)
        for row in SQ:
            header.append(row)
    header.append('@PG\tID:proBamPy\tVN:1.0\tCL:'+str(command_line))
    header.append('@GA\tAS:'+str(database)+'\tVN:'+str(database_v))
    # get comments and append comments to file
    if comments!=[]:
        for comment in comments:
            comment=str(comment).rstrip()
            comment=re.sub(r'(^[ \t]+|[ \t]+(?=:))', '', comment, flags=re.M)
            header.append('@CO\t'+str(comment))
    comments=extract_comments(psm_file)
    if comments!=[]:
        for comment in comments:
            comment=str(comment).rstrip()
            comment=re.sub(r'(^[ \t]+|[ \t]+(?=:))', '', comment, flags=re.M)
            header.append('@CO\t'+str(comment))
    for row in header:
        file.write(row+'\n')


#
# Function to convert SAM to BAM
#
def sam_2_bam(directory,name):
    '''
    :param directory:
    :param name: file name
    '''
    print("Converting SAM to BAM")
    infile = pysam.AlignmentFile((directory+name+'.sam'), "r")
    outfile = pysam.AlignmentFile((directory+name+'.bam'), "wb", template=infile)
    for s in infile:
        outfile.write(s)
    infile.close()
    outfile.close()
    # create EOF
    bam=open((directory+name+'.bam'),'ab')
    bam.write("\x1f\x8b\x08\x04\x00\x00\x00\x00\x00\xff\x06\x00BC" + \
                "\x02\x00\x1b\x00\x03\x00\x00\x00\x00\x00\x00\x00\x00\x00")
    bam.close()

    # Pysam v 0.8.4.:
    pysam.sort((directory + name + '.bam'), (directory + name + '.sorted'))

    # For new pysam version, has error for bigger files
    # pysam.sort("-o",(directory+name+'.sorted.bam'),(directory+name+'.bam'))
    pysam.index(directory+name+'.sorted.bam')
#
# function to calculate and adjust NH for every peptide
# Also depending on mode, can create peptide based proBAM
#

def compute_NH_XL(directory,name,include_unmapped,mode):
    '''
    :param directory: directory of the proBAM file
    :param name: proBAM file name
    :param include_unmapped: 'Y' or 'N', whether to include unmapped PSMs
    :param mode: proBAM mode (psm,peptide,peptide-mod...)
    '''
    if mode in ['proBAM_peptide', 'proBAM_peptide_mod']:
        print("Create peptide-based proBAM")
        sam_file = open(directory + name + '.sam', 'r')
        original_file = sam_file.read()
        nh_hash = {}
        score_hash = {}
        peptide_hash={}
        for line in original_file.split("\n"):
            if len(line) < 1:
                continue
            elif line[0] == "@":
                continue
            else:
                if line.split("\t")[5] == "*":
                    continue
                else:
                    line=line.split("\t")
                    if mode=='proBAM_peptide':
                        key=line[23]+'_'+line[2]+'_'+line[3]+'_'+line[5]
                        id=line[23].split(':')[2]
                    elif mode=='proBAM_peptide_mod':
                        if line[20]=="XM:Z:*":
                            id=line[23].split(':')[2]
                        else:
                            id=line[23].split(':')[2]+","+line[20].split('XM:Z:')[1].replace(";",",")
                        key=id+'_'+line[2]+'_'+line[3]+'_'+line[5]

                    if id in nh_hash:
                        if create_id_from_list([line[2], line[3], line[5]]) in \
                                nh_hash[id]:
                            continue
                        else:
                            nh_hash[id].append(
                                create_id_from_list([line[2], line[3],
                                                     line[5]]))
                    else:
                        nh_hash[id] = [(create_id_from_list([line[2], line[3],
                                                                            line[5]]))]
                    if key not in peptide_hash:
                        peptide_hash[key]=[id,line[1],line[2],line[3],line[4],line[5],line[6],line[7],line[8],line[9],
                                           line[10],line[11],line[12],'XB:f:*','XC:i:*',line[15],line[16],line[17],
                                           'XI:f:*','XL:i:*','XM:Z:*',line[21],line[22],line[23],line[24],line[25],
                                           line[26],line[27],line[28],line[29],line[30],line[31],line[32]]
                    if id not in score_hash:
                        score_hash[id]=1
                    else:
                        if line[26]!='XS:f:*':
                            try:
                                if float(line[26].split('XS:f:')[1])>float(peptide_hash[key][26].split('XS:f:')[1]):
                                    peptide_hash[key][26]=line[26]
                            except:pass
                        if line[24] != 'XQ:f:*':
                            try:
                                if float(line[24].split('XQ:f:')[1]) > float(peptide_hash[key][24].split('XQ:f:')[1]):
                                    peptide_hash[key][24] = line[24]
                            except:
                                pass

        sam_file.close()
        sam_file = open(directory + name + '.sam', 'w')
        for line in original_file.split("\n"):
            if len(line) < 1:
                continue
            elif line[0] == "@":
                sam_file.write(line)
                sam_file.write("\n")
        for key in peptide_hash:
            if include_unmapped=='N' and peptide_hash[key][5]=='*':
                continue
            else:
                peptide_hash[key][11] = "NH:i:"+ str(len(nh_hash[peptide_hash[key][0]]))
                sam_file.write(peptide_hash[key][0])
                for col in peptide_hash[key][1:]:
                    sam_file.write('\t'+col)
                sam_file.write("\n")
    else:
        sam_file=open(directory+name+'.sam','r')
        original_file = sam_file.read()
        nh_hash={}
        xl_hash={}
        for line in original_file.split("\n"):
            if len(line)<1:
                continue
            elif line[0]=="@":
                continue
            else:
                if line.split("\t")[5]=="*":
                    continue
                else:
                    if line.split("\t")[0] in xl_hash:
                        if line.split("\t")[23] not in xl_hash[line.split("\t")[0]]:
                            xl_hash[line.split("\t")[0]].append(line.split("\t")[23])
                    else:
                        xl_hash[line.split("\t")[0]]=[]
                        xl_hash[line.split("\t")[0]].append(line.split("\t")[23])
                    if nh_key_line(line) in nh_hash:
                        if create_id_from_list([line.split('\t')[2],line.split('\t')[3],line.split('\t')[5]]) in \
                                nh_hash[nh_key_line(line)]:
                            continue
                        else:
                            nh_hash[nh_key_line(line)].append(create_id_from_list([line.split('\t')[2],line.split('\t')[3],
                                                                                   line.split('\t')[5]]))
                    else:
                        nh_hash[nh_key_line(line)]=[(create_id_from_list([line.split('\t')[2], line.split('\t')[3],
                                                                          line.split('\t')[5]]))]
        sam_file.close()
        sam_file=open(directory+name+'.sam','w')
        for line in original_file.split("\n"):
            if len(line)<1:
                continue
            elif line[0]=="@":
                sam_file.write(line)
            elif line.split("\t")[5]=="*":
                if include_unmapped!='N':
                    sam_file.write(line)
                else:
                    continue
            else:
                line=line.replace("XL:i:*","XL:i:"+str(len(xl_hash[line.split("\t")[0]])))
                line=line.replace("NH:i:*","NH:i:"+str(len(nh_hash[nh_key_line(line)])))
                sam_file.write(line)
            sam_file.write("\n")

####################
### MAIN PROGRAM ###
####################

if __name__ == '__main__':
    if len(sys.argv) > 1:
        parser = get_parser()  # start command line argument parser
        args = parser.parse_args()  # parse command line arguments
        globals().update(args.__dict__)
        start_time = time.time()
        # start timing function
        get_input_variables()
        database_v=int(database_v)
        allowed_mismatches=int(allowed_mismatches)
        directory=directory+'/'

        # hash PSM_DATA and define variables
        psm_hash=proBAM_input.get_PSM_hash(psm_file,decoy_annotation)

        parse_results=proBAM_IDparser.parseID(psm_hash,species,database,decoy_annotation,database_v,three_frame_translation
                                              ,pre_picked_annotation)
        annotation=parse_results[1]
        psm_hash=parse_results[0]
        transcript_hash=annotation[0]
        exon_hash=annotation[1]
        id_map=parse_results[2]

        # convert to SAM
        if conversion_mode!='proBED':
            file = open_sam_file(directory, name)
            create_SAM_header(file,version, database, sorting_order, database_v, species, command_line, psm_file, comments)
            PSM2SAM(psm_hash,transcript_hash,exon_hash,decoy_annotation,allowed_mismatches,file,rm_duplicates,
                    three_frame_translation,psm_file,id_map,None)
            compute_NH_XL(directory, name, include_unmapped,conversion_mode)
            sam_2_bam(directory, name)
        # convert to BED
        else:
            file = proBAM_proBED.open_bed_file(directory, name)
            proBAM_proBED.create_BED_header(file, database, database_v, command_line, psm_file, comments)
            proBAM_proBED.PSM2BED(psm_hash,transcript_hash,exon_hash,decoy_annotation,allowed_mismatches,file,
                                  rm_duplicates,three_frame_translation,id_map,None,database_v,species)


        print("proBAM conversion succesful")
        print("%f seconds" % (time.time() - start_time))         # output script run time
    else:
        # start GUI
        import proBAM_GUI
        proBAM_GUI._get_global_arguments_()
        proBAM_GUI.GUI()
