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
from itertools import imap
import operator
import re
import proBAM_pepxml
import proBAM_mzTab
import proBAM_mzid
from cogent.core.genetic_code import DEFAULT as standard_code
#
# reverse complements a DNA string
#
def reverse_complement(DNA):
    '''
    :param DNA: DNA string
    :return: reverse complement of the DNA string
    '''
    temp=DNA[::-1]
    DNA=""
    for i in temp:
        if i=="A":
            DNA+='T'
        elif i=="T":
            DNA+='A'
        elif i=="G":
            DNA+='C'
        else:
            DNA+='G'
    return DNA

#
# Function to calculate distance between 2 strings
#
def hamming(str1, str2):
    '''
    :param str1: first string
    :param str2: second string
    :return: the hamming distance
    '''
    assert len(str1) == len(str2)
    #ne = str.__ne__  ## this is surprisingly slow
    ne = operator.ne
    return sum(imap(ne, str1, str2))

#
# Function to compute the CIGAR string of a protein
#
def compute_cigar(gen_pos,exons,strand,peptide):
    '''
    :param gen_pos: genomic start position
    :param exons: exon dictionairy
    :param strand: transcript strand
    :param peptide: peptide sequence
    :return: CIGAR string
    '''
    temp_exons=exons
    gen_pos=int(gen_pos)
    hit=0
    cigar=''
    frame_pos=''
    pointer=0
    debug=0
    # debug section for wrong CIGAR
    #if peptide=="SCCSCCPVGCAK":
    #    debug=1
    length= int(len(peptide)*3)
    temp_exons=sorted(temp_exons,key=lambda x: int(x[2]))
    if strand=='-1':
        adjusted_temp_exons=list(reversed(temp_exons))
    else:
        adjusted_temp_exons=temp_exons
    for exon in adjusted_temp_exons:
        exon=[int(numeric_string) for numeric_string in exon]
        if hit ==0:
            if gen_pos > exon[1]:
                continue
            else:
                if (gen_pos+length)<=exon[1]+1:
                    cigar=str(length)+'M'
                    frame_pos+=str(gen_pos%3)
                    break
                else:
                    pointer=(exon[1]-gen_pos+1)
                    cigar=cigar+(str(pointer)+'M')
                    length=length-pointer
                    prev_exon_end=exon[1]
                    prev_rank=exon[2]
                    hit=1
                    frame_pos += str(gen_pos % 3)+","
                    continue
        if hit!=0:
            cigar+=(str((exon[0]-prev_exon_end-1))+'N')
            if '-' in cigar:
                print strand
                print peptide
                print cigar
                print adjusted_temp_exons
                print gen_pos, length
            if (exon[0]+(length)-1)>exon[1]+1:
                cigar=cigar+(str((exon[1]-exon[0]+1))+'M')
                length=length-(exon[1]-exon[0]+1)
                frame_pos+=str(exon[0]%3)+','
                prev_exon_end = exon[1]
                continue
            else:
                frame_pos+=str(exon[0]%3)
                cigar=cigar+(str(length)+'M')
                break
    #debug section fro wrong CIGAR
    '''
    if len(cigar.split('M'))>3:
        print strand
        print peptide
        print cigar
        print temp_exons
        print gen_pos,length
    '''
    return [cigar,frame_pos]
#
# Function returns genomic position of the leftmost base
#
#
def calculate_genome_position(phit,strand,offset,start_exon_rank,peptide,exons,chr,three_frame_translation,shift):
    '''
    :param phit: location of peptide hit
    :param strand: transcript strand
    :param offset: CDS offset
    :param start_exon_rank: rank of first exon containing the CDS
    :param peptide: peptide sequence
    :param exons: exon hash
    :param chr: transcript chr
    :return: genomic start position
    '''
    if three_frame_translation!='Y':
        temp_offset=offset-shift
    temp_exons=exons
    if three_frame_translation=='Y':
        if strand=='1':
            tr_pos=(phit)
            pointer=0
            #start_exon_rank=0
        else:
            tr_pos=(phit)+(len(peptide)*3)-1
            #start_exon_rank=len(temp_exons)+1
            pointer=0
        gen_pos=0
        temp_offset=0
        start_exon_rank=1

    else:
        if strand=='1':
            tr_pos=(phit*3)
        else:
            tr_pos=(phit*3)+(len(peptide)*3)-1
        pointer=1
        gen_pos=0

    # iterate over temp_exons till peptide position is encoutered
    temp_exons=sorted(temp_exons,key=lambda x: int(x[2]))
    for exon in temp_exons:
        exon=[int(numeric_string) for numeric_string in exon]
        if exon[2]<start_exon_rank:
            continue
        elif exon[2]==start_exon_rank:
            if strand=='1':
                exon[0]=exon[0]+temp_offset
            else:
                exon[1]=exon[1]-temp_offset
        if tr_pos>((exon[1]-exon[0])+pointer):
            pointer=pointer+(exon[1]-exon[0]+1)
            continue
        else:
            remember_exon=exon
            if strand=='1':
                gen_pos=exon[0]+((tr_pos-pointer))
            else:
                gen_pos=exon[1]-((tr_pos-pointer))
            break

    # below print tests are generated to get the left most basepare
    '''
    print '*************************************************'
    print 'test:',start_exon_rank,phit,strand
    print pointer,tr_pos, temp_offset
    #print remember_exon
    print gen_pos,peptide,strand,chr
    print '*************************************************'
    print ''
    import sys
    sys.exit()
    '''

    return [str(gen_pos),temp_exons]

#
#Function that maps peptide on the corresponding protein, alowing mismatches ( as specified)
# and returns the first left base position mapped after translating the transcript sequence in 3-frames
#
def map_peptide_to_protein_3frame(peptide_seq,transcript_seq,allowed_mismatches,strand):
    '''
    :param peptide_seq: peptide sequence (string)
    :param transcript_seq: transcript sequence (string)
    :param allowed_mismatches: number of allowed mismatches
    :param strand: chromosome strand
    :return: number of hits of peptide on protein
    '''
    size_adjust=-1    # adjust size of transcript for starting at +1/+2frame
    hits=[]
    pre_post_aa=['','']
    pep_length=len(peptide_seq)
    frame=[0]*3
    frame[0]=standard_code.translate(transcript_seq)
    frame[1]=standard_code.translate(transcript_seq[1:])
    frame[2]=standard_code.translate(transcript_seq[2:])
    for f in frame:
        size_adjust+=1
        for i in range(0,(len(f)-pep_length)):
            if hamming(peptide_seq,f[i:pep_length+i]) <= allowed_mismatches:
                adjusted_hit_pos=(i*3)+size_adjust
                hits.append([adjusted_hit_pos,hamming(peptide_seq,f[i:pep_length+i])])

                # compute 2 preceding AA
                if (i - 1) == 0:
                    pre_post_aa[0] = f[(i - 1)]
                elif (i - 2) >= 0:
                    pre_post_aa[0] = f[(i - 2):i]
                else:
                    pre_post_aa[0] = "*"


                # compute 2 folowwing AA
                if (i + 1) == (len(f) - 1):
                    pre_post_aa[1] = f[pep_length + i]
                elif (i + 2) <= (len(f) - 1):
                    pre_post_aa[1] = f[(pep_length + i):(pep_length + i + 2)]
                else:
                    pre_post_aa[1] = "*"
    return [hits,pre_post_aa]
#
# Function that maps peptide on the corresponding protein, alowing mismatches ( as specified)
# and returns the first left base position mapped
#
def map_peptide_to_protein(peptide_seq,protein_seq,allowed_mismatches):
    '''
    :param peptide_seq: peptide sequence
    :param protein_seq: protein sequence
    :param allowed_mismatches: allowed mismatches
    :return: location of mapped peptide
    '''
    hits=[]
    pep_length=len(peptide_seq)
    prot_length=len(protein_seq)
    pre_post_aa=['','']
    for i in range(0,prot_length-pep_length):
        if hamming(peptide_seq,protein_seq[i:pep_length+i]) <= allowed_mismatches:
            hits.append([i,hamming(peptide_seq,protein_seq[i:pep_length+i])])
            # compute 2 preceding AA
            if (i - 1) == 0:
                pre_post_aa[0] = protein_seq[(i - 1)]
            elif (i - 2) >= 0:
                pre_post_aa[0] = protein_seq[(i - 2):i]
            else:
                pre_post_aa[0] = "*"
            # compute 2 folowwing AA
            if (i + 1) == (len(protein_seq) - 1):
                pre_post_aa[1] = protein_seq[pep_length + i]
            elif (i + 2) <= (len(protein_seq) - 1):
                pre_post_aa[1] = protein_seq[(pep_length + i):(pep_length + i + 2)]
            else:
                pre_post_aa[1]="*"
    return [hits,pre_post_aa]

#
# Function to calculate FLAG of PSM
#
def calculate_FLAG(strand,rank,unmapped):
    '''
    :param strand: transcript strand
    :param rank: psm rank
    :param decoy: boolean decoy
    :return: FLAG
    '''
    strand=int(strand)
    FLAG=0
    if strand==-1:
        FLAG=FLAG+16
    if str(rank)!=str(1) and rank!='*' and rank!=[] and rank!="":
        FLAG=FLAG+256
    if unmapped==1:
        FLAG=FLAG+4
    return FLAG

#
# Function to analyze modification and create XA-string
#
def create_XM(modifications):
    '''
    :param modifications: psm modification
    :return: modification in proBAM format
    '''
    XM=''
    if modifications==[]:
        return'*'
    else:
        for mod in modifications:
            if XM!='':
                XM+=';'
            XM+=str(mod['position'])
            XM+='-'
            XM+=str(mod['mass'])
    return XM

#
# Function to determine peptide annotation and create XA string
#
def create_XA(mismatch):
    '''
    :param mismatch: number of allowed mismatches
    :return: XA
    '''
    XA='XA:i:'
    if mismatch==0:
        XA+='0'
    else:
        XA+='-1'
    return XA
#
# Function to analyze peptide type and create XG string
#
def create_XG(hamming_distance):
    '''
    :param hamming_distance: hamming distance to use
    :return: the hamming distance
    '''
    XG='XG:Z:'
    if hamming_distance==0:
        XG+='N'
    else:
        XG+='V'
    return XG

#
# get enzyme
#
def get_enzyme(psm_file):
    '''
    :param psm_file: psm file
    :return: enzyme used in the psm file
    '''
    #catch mzid
    if re.match('^.*\.(mzid)$',psm_file.lower())!=None:
        return proBAM_mzid.get_enzyme_mzid(psm_file)
    # catch pepxml file format and parse
    elif re.match('^.*\.(pepxml|pep.xml|xml)$',psm_file.lower())!=None:
        return proBAM_pepxml.get_enzyme_pepxml(psm_file)
    # catch mztab file format and parse
    elif re.match('^.*\.(mztab)$',psm_file.lower())!=None:
        return proBAM_mzTab.get_enzyme_mztab(psm_file)

#
# get enzyme specificity
#
def get_enzyme_specificity(psm_file):
    '''
    :param psm_file: psm file
    :return: enzyme specificity extracted from psm_file
    '''
    # catch mzid
    if re.match('^.*\.(mzid)$', psm_file.lower()) != None:
        return proBAM_mzid.get_enzyme_specificity_mzid(psm_file)
    # catch pepxml file format and parse
    elif re.match('^.*\.(pepxml|pep.xml|xml)$', psm_file.lower()) != None:
        return proBAM_pepxml.get_enzyme_specificity_pepxml(psm_file)
    # catch mztab file format and parse
    elif re.match('^.*\.(mztab)$', psm_file.lower()) != None:
        return proBAM_mzTab.get_enzyme_specificity_mztab(psm_file)
#
# translate RNA sequence
#
def translate_seq(sequence,strand):
    '''
    :param sequence: DNA sequence
    :param strand: strand of DNA sequence
    :return: translated DNA sequence
    '''
    if str(strand)=="1":
        return standard_code.translate(sequence)
    else:
        return standard_code.translate(reverse_complement(sequence))
#
# Extract comments containing experiment information for PSM file
#
def extract_comments(psm_file):
    '''
    :param psm_file: psm file
    :return: returns an array of extracted comments from the original psm file
    '''
    if re.match('^.*\.(mzid)$', psm_file.lower()) != None:
        return proBAM_mzid.extract_comments_from_mzid(psm_file)
    # catch pepxml file format and parse
    elif re.match('^.*\.(pepxml|pep.xml|xml)$', psm_file.lower()) != None:
        return proBAM_pepxml.extract_comments_from_pepxml(psm_file)
    # catch mztab file format and parse
    elif re.match('^.*\.(mztab)$', psm_file.lower()) != None:
        return proBAM_mzTab.extract_comments_from_mztab(psm_file)
#
# Give a list, create a identifiers where elements from the list are concatenated by '_'
#
def create_id_from_list(list):
    '''
    :param list: list of elements (strings)
    :return: identifiers
    '''
    id=""
    for i in list:
        if id == "":
            id+=str(i)
        else:
            id+="_"+str(i)
    return id
#
# for a line creates a unique genomic location key for this peptide
#
def nh_key_line(line):
    '''
    :param line:a line from a SAM file
    :return: key (string)
    '''
    line=line.split("\t")
    key=line[11]+"_"+line[0]
    return key
