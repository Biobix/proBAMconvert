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

__author__ = 'vladie'

from pyteomics import mzid
from fnmatch import fnmatch
import re

#
# Parse mzid file into dictionairy
#

def get_PSM_mzid(psm_file,validated_only):
    '''
    :param psm_file: mzid file
    :return: dictionairy of parsed mzid file, suitable for proBAMconvert
    '''

    with mzid.read(psm_file) as PSM:
        psm_hash=[]
        accession_hash=_get_accessions_(psm_file)
        mod_hash=_get_modification_(psm_file)
        sequence_hash=_get_peptide_sequence_hash(psm_file)
        spectraData_ref={}
        c=0
        for row in PSM:
            if 'spectraData_ref' in row:
                if row['spectraData_ref'] not in spectraData_ref:
                    spectraData_ref[row['spectraData_ref']]=c
                    c+=1
                row['spectrumID']='ms_run['+str(spectraData_ref[row['spectraData_ref']])+']:'+row['spectrumID']
            temp_hash={"assumed_charge":row['SpectrumIdentificationItem'][0]['chargeState'],"spectrum":row['spectrumID'],"search_hit":[]}
            for psm in row["SpectrumIdentificationItem"]:
                if psm['passThreshold']==True or validated_only=="N":
                    proteins=[]
                    massdiff=_cal_massdiff_(psm['experimentalMassToCharge'],psm['calculatedMassToCharge'])
                    for protein in psm["PeptideEvidenceRef"]:
                        proteins.append({"protein":accession_hash[protein['peptideEvidence_ref']]})
                    temp_hash['search_hit'].append({"hit_rank":psm['rank'],"modifications":mod_hash[psm['peptide_ref']],
                                                    "calc_neutral_pep_mass": psm['experimentalMassToCharge'],
                                                    "precursor_neutral_mass": psm['calculatedMassToCharge'],
                                                    "peptide":sequence_hash[psm['peptide_ref']],
                                                    "search_score":{"score":_get_score_(psm),"evalue":_get_evalue_(psm)},
                                                    "proteins":proteins,"num_missed_cleavages":"0",
                                                    "massdiff":massdiff})
            psm_hash.append(temp_hash)
    PSM.close()
    del mod_hash
    del sequence_hash
    del accession_hash
    return psm_hash
#
# Depreciated
#
def _filter_peptide_ref_(seq):
    '''
    :param seq: peptide sequence
    :return: filtered peptides sequence (just the sequence no modication or pre/postfixes
    '''
    if '_' in seq:
        seq.split("_")[0]
    return re.sub("[^a-zA-Z]+", "", seq)
#
# Retrieve modifications from the psm file
#
def _get_modification_(psm_file):
    '''
    :param psm_file: psm file
    :return: dictionairy of peptide modification for every peptide sequence
    '''
    mod_hash = {}
    with open(psm_file, 'r') as f:
        hit=0
        for line in f:

            if "<Peptide " in line:
                mod=[]
                uni_mod=[]
                temp_line = line.replace('><', ' ')
                temp_line = temp_line.replace("/>", '')
                temp_line = temp_line.replace(">", '')
                temp_line = temp_line.replace('\n', '')
                temp_line = temp_line.replace('\r','')
                temp_line = temp_line.split(" ")
                for tag in temp_line:
                    tag = tag.split("=")
                    if tag[0] == "id":
                        id = tag[1].split(">")[0].replace("\"", "")
                hit=1

            if hit==1:
                loc=0
                mass=0
                temp_line = line.replace("/>", '')
                temp_line = temp_line.replace(">", '')
                temp_line = temp_line.replace('\n', '')
                temp_line = temp_line.split(" ")
                for tag in temp_line:
                    tag = tag.split("=")
                    if tag[0]=="location":
                        loc=tag[1].split(">")[0].replace("\"","")
                        if mass!=0:
                            mod.append({"position": loc, "mass": mass})
                    if tag[0]=="monoisotopicMassDelta":
                        mass=tag[1].split(">")[0].replace("\"","")
                        if loc!=0:
                            mod.append({"position":loc,"mass":mass})
                    if tag[0]=="accession" and ("MOD" in tag[1]):
                        tag[1]=tag[1].split(">")[0].replace("\"","")
                        mass=tag[1]
                        uni_mod.append({"position":loc,"mass":mass})

            if "/Peptide>" in line:
                if uni_mod==[]:
                    mod_hash[id]=mod
                else:
                    mod_hash[id]=uni_mod
                hit=0
    return mod_hash
#
# Retrieve peptide sequences from the psm file
#
def _get_peptide_sequence_hash(psm_file):
    '''
    :param psm_file: psm file
    :return: dictionairy of peptide sequences
    '''
    sequence_hash = {}
    with open(psm_file, 'r') as f:
        hit = 0
        for line in f:

            if "<Peptide " in line:
                mod = []
                uni_mod = []
                temp_line = line.replace('><', ' ')
                temp_line = temp_line.replace("/>", '')
                temp_line = temp_line.replace(">", '')
                temp_line = temp_line.replace('\n', '')
                temp_line = temp_line.replace('\r', '')
                temp_line = temp_line.split(" ")
                for tag in temp_line:
                    tag = tag.split("=")
                    if tag[0] == "id":
                        id = tag[1].split(">")[0].replace("\"", "")
                hit = 1

            if hit == 1:
                loc = 0
                mass = 0
                if '<PeptideSequence>' in line:
                    temp_line=line.split('<PeptideSequence>')[1]
                    sequence=temp_line.split('</PeptideSequence>')[0]
            if "/Peptide>" in line:
                sequence_hash[id] = sequence
                hit = 0
        return sequence_hash
#
# Extract accession from psm file (contains the protein identifiers)
#
def _get_accessions_(psm_file):
    '''
    :param psm_file: psm file
    :return: dictionaries with accessions
    '''
    accession_hash={}
    peptide_evidence_hash={}
    dbsequence_hash={}
    with open(psm_file,'r') as f:
        for line in f:
            if "<PeptideEvidence" in line:
                if "id" in line and "dBSequence_ref" in line:
                    line=line.replace("/>",'')
                    line=line.replace(">",'')
                    line=line.replace('\n','')
                    line=line.split("</")[0]
                    line=line.split(" ")
                    accession=""
                    for tag in line:
                        tag=tag.split("=")
                        if tag[0]=="dBSequence_ref":
                            dbseq_ref=tag[1].replace("\"","")
                        elif tag[0]=="id":
                            id=tag[1].replace("\"","")
                        elif tag[0] == "accession":
                            dbseq_ref = tag[1].replace("\"", "")
                    if dbseq_ref not in peptide_evidence_hash:
                        peptide_evidence_hash[dbseq_ref]=[id]
                    else:
                        peptide_evidence_hash[dbseq_ref].append(id)
            elif "<DBSequence" in line:
                if "id" in line and "accession" in line:
                    line = line.replace("/>", '')
                    line = line.replace(">", '')
                    line = line.replace('\n', '')
                    line = line.split("</")[0]
                    line = line.split(" ")
                    for tag in line:
                        tag = tag.split("=")
                        if tag[0] == "accession":
                            accession = tag[1].replace("\"", "")
                        elif tag[0] == "id":
                            id = tag[1].replace("\"", "")
                    dbsequence_hash[id] = accession
    for key in dbsequence_hash:
        if key in peptide_evidence_hash:
            for id in peptide_evidence_hash[key]:
                accession_hash[id]=dbsequence_hash[key]
    for key in peptide_evidence_hash:
        for id in peptide_evidence_hash[key]:
            if id not in accession_hash:
                accession_hash[id]=key
    return accession_hash
#
# Calculate the mass difference given the calculated and experimental mass
#
def _cal_massdiff_(calc_mass,exp_mass):
    '''
    :param calc_mass: calculated mass
    :param exp_mass: experimental mass
    :return: mass difference
    '''
    mass_diff=exp_mass-calc_mass
    return str(mass_diff)+';'+str(exp_mass)+';'+str(calc_mass)
#
# Extract the score for a psm
#
def _get_score_(psm):
    """
    :param psm: peptide to spectrum match dictionairy
    :return: XCorr score (if not available XCorr = 0)
    """
    hit=0
    hit_key=''
    score={}
    #print psm.keys()
    for key in psm:
        if "score" in key.lower():
            hit=1
            hit_key=key
    if hit==1:
        return psm[hit_key]
    else:
        return "-1"
#
# Extract the e-value for a psm
#
def _get_evalue_(psm):
    '''
    :param psm: psm
    :return: the e-value if any
    '''
    hit=0
    hit_key=''
    score={}
    #print psm.keys()
    for key in psm:
        if "xcorr" in key.lower() or 'expectation' in key.lower() or 'confidence' in key.lower() \
        or "e_value" in key.lower().replace("-","_") or 'evalue' in key.lower() or 'fdr' in key.lower():
            hit=1
            hit_key=key
    if hit==1:
        return psm[hit_key]
    else:
        return "-1"
#
# Create the modified peptide sequence
#
def _get_mod_peptide_sequence_(sequence,modification):
    '''
    :param sequence: peptide sequence with optional tags
    :return: peptide sequence without tags
    '''
    mod_peptide_sequence=""
    for position in range(0,len(sequence)):
        mod_peptide_sequence+=sequence[position]
        for mod in modification:
            if int(mod['location'])==position:
                mod_peptide_sequence+='['+str(mod['monoisotopicMassDelta'])+']'

    return mod_peptide_sequence
#
# get modification for a psm and rearrange into right format
#
def _get_peptide_modifications_(modifications):
    '''
    :param sequence: peptide sequence
    :return: modification parsed out peptide sequence
    '''
    modification=[]
    for mod in modifications:
        pos=mod['location']
        value=mod['monoisotopicMassDelta']
        modification.append({"position":pos,"mass":value})
    return modification
#
# Function to extract enzyme from psm_file
#
def get_enzyme_mzid(psm_file):
    '''
    :param psm_file: psm file
    :return: enzyme if able to extract from psm_file
    '''
    f=open(psm_file,'r')
    lines = f.read()
    start = lines.find('<EnzymeName>')
    stop = lines.find('</EnzymeName>')
    line = lines[start:stop]
    line=line.replace('-','_')
    line=line.lower()
    if fnmatch(line,'*trypsin?p*'):
        return 2
    elif fnmatch(line,'*trypsin*'):
        return 1
    elif fnmatch(line, '*lys?c*'):
        return 3
    elif fnmatch(line, '*lys?n*'):
        return 4
    elif fnmatch(line, '*arg?c*'):
        return 5
    elif fnmatch(line, '*asp?n*'):
        return 6
    elif fnmatch(line, '*cnbr*'):
        return 7
    elif fnmatch(line, '*glu?c*'):
        return 8
    elif fnmatch(line, '*pepsina*'):
        return 9
    elif fnmatch(line, '*chymotrypsin*'):
        return 10
    elif fnmatch(line, '*noenzyme*'):
        return 0
    else:
        return "-1"
#
# Attempt to extract the enzyme specificity from the psm_file
#
def get_enzyme_specificity_mzid(psm_file):
    '''
    :param psm_file: psm file
    :return: enzyme specificity
    '''
    f=open(psm_file,'r')
    lines = f.read()
    start = lines.find('<EnzymeName>')
    stop = lines.find('</EnzymeName>')
    line = lines[start:stop]
    line=line.replace('-','_')
    line=line.lower()
    if fnmatch(line,'*semiSpecific=1*') or fnmatch(line,'*semiSpecific=true*'):
        return 1
    elif fnmatch(line,'*semiSpecific=0*') or fnmatch(line,'*semiSpecific=false*'):
        return 2
    else:
        return 3
#
# Ectract comments from the psm_file
#
def extract_comments_from_mzid(psm_file):
    '''
    :param psm_file: psm file
    :return: comments
    '''
    f=open(psm_file,'r')
    comments=[]
    lines = f.read()

    start = lines.find('<AnalysisSoftwareList')
    stop = lines.find('</AnalysisSoftwareList')
    line = lines[start:stop]
    for comment in line.split("\n"):
        comments.append(comment)

    start = lines.find('<AnalysisProtocolCollection')
    stop = lines.find('</AnalysisProtocolCollection')
    line = lines[start:stop]
    for comment in line.split("\n"):
        comments.append(comment)

    start = lines.find('<AnalysisCollection')
    stop = lines.find('</AnalysisCollection')
    line = lines[start:stop]
    for comment in line.split("\n"):
        comments.append(comment)
    return comments