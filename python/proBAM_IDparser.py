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

__author__="vladie"
import proBAM_ENSEMBL
import proBAM_biomart
import re
import sys
import urllib,urllib2
from bioservices import UniProt

#
# Find accession encryption pattern and database
#

def parseID(psm_hash,species,database,decoy_annotation,database_v,three_frame_translation,pre_picked_annotation):
    '''
    :param psm_hash: psm hash (see: proBAM_input.get_PSM_hash())
    :param species: species
    :param database: annotation database name
    :param decoy_annotation: decoy annotation list
    :param database_v: annotation database version
    :return: [PSM_hash_with_updated_identifiers,[transcript_annotation_hash,exon_hash]]
    '''
    # create list with IDs, for faster access
    if pre_picked_annotation=='all':
        return parse_all_IDs(psm_hash,species,database,decoy_annotation,database_v,three_frame_translation,
                             pre_picked_annotation)
    else:
        PID_list=[]
        for key in psm_hash:
                if 'search_hit' in key:
                    for psm in key['search_hit']:
                        for i in range (0,len(psm['proteins'])):
                           PID_list.append(psm['proteins'][i]['protein'])
        found=0
        count=0
        # set maximum amount of identifiers to be checked before giving up identifier identification
        if len(PID_list)<1000:
            max=len(PID_list)
        else:
            max=1000
        # parse | and _ out of identifier and put in a list of terms
        decryption=[]
        while found==0 and count<max:
            find_result = _find_annotation_(PID_list[count],species, pre_picked_annotation)
            if find_result[0]==1:
                found=1
                break
            count+=1
        # print PID_list[0:100]
        # raise error max amount of identifiers checked and nog identifiers identified
        if found==0:
            raise ValueError('Protein ID annotations not recognized\n'\
                'currently supported Protein ID annotations:\n'\
                'ENSEMBL protein identifiers\n ENSEMBL transcript identifiers\n UNIPROT/SWISSPROT identifiers '
                 '\n UNIPROT ENTRIES \n REFSEQ identifiers')
            sys.exit()

        # if identifiers identified, update PSM_hash to the right identifiers
        elif found==1:
            protein_ID=[]
            for key in psm_hash:
                if 'search_hit' in key:
                    for psm in key['search_hit']:
                        for i in range (0,len(psm['proteins'])):
                            psm['proteins'][i]['protein']=_update_protein_accession_(psm['proteins'][i]['protein'],
                                                                                   decoy_annotation,find_result[2])
                            hit=0
                            if 'DECOY_' in psm['proteins'][i]['protein'].upper():
                                hit=1
                            if hit==0:
                                protein_ID.append(psm['proteins'][i]['protein'])

            # if database == ENSEMBL use correct function to retrieve anntotation information
            if database=="ENSEMBL":
                if find_result[1]=="ENSEMBL_TR":
                    id_map = _id_map_('ENSEMBL', 'ENSEMBL', protein_ID, psm_hash, species, decoy_annotation,database_v)
                    annotation= proBAM_ENSEMBL.prepareAnnotationENSEMBL(_get_transcript_ids_from_map(id_map),'transcript',
                                                                        database_v,species,three_frame_translation)
                elif find_result[1]=="ENSEMBL_PR":
                    id_map = _id_map_('ENSEMBL', 'ENSEMBL', protein_ID, psm_hash, species, decoy_annotation,database_v)
                    annotation= proBAM_ENSEMBL.prepareAnnotationENSEMBL(_get_transcript_ids_from_map(id_map),'protein',
                                                                        database_v,species,three_frame_translation)
                elif find_result[1]=="UNIPROT":
                    id_map = _id_map_('UNIPROT','ENSEMBL',protein_ID,psm_hash,species,decoy_annotation,database_v)
                    annotation= proBAM_ENSEMBL.prepareAnnotationENSEMBL(_get_transcript_ids_from_map(id_map),'transcript',
                                                                        database_v,species,three_frame_translation)
                elif find_result[1] == "UNIPROT_ENTRY":
                    id_map = _id_map_('UNIPROT_ENTRY', 'ENSEMBL', protein_ID, psm_hash, species, decoy_annotation,
                                          database_v)
                    annotation = proBAM_ENSEMBL.prepareAnnotationENSEMBL(_get_transcript_ids_from_map(id_map), 'transcript',
                                                                         database_v, species, three_frame_translation, )
                elif find_result[1]=="REFSEQ":
                    id_map = _id_map_('REFSEQ', 'ENSEMBL', protein_ID, psm_hash, species, decoy_annotation,
                                      database_v)
                    annotation = proBAM_ENSEMBL.prepareAnnotationENSEMBL(_get_transcript_ids_from_map(id_map),
                                                                         'transcript',
                                                                         database_v, species, three_frame_translation, )
                return [psm_hash,annotation,id_map]
            # raise error if database unsupported
            else:
                raise ValueError('Currently supported annotation databases: \n ENSEMBL')
#
# Parse all ID's (multiple ID annotations)
#
def parse_all_IDs(psm_hash,species,database,decoy_annotation,database_v,three_frame_translation,pre_picked_annotation):
    '''
    :param psm_hash: psm dictionairy
    :param species: species
    :param decoy_annotation: which decoy annotation is used
    :param database_v: database version
    :param three_frame_translation: whether transcript should be 3-frame translated
    :return: array with updated psm dictionairy, and mapped annotations and identifiers to Ensembl
    '''
    ensembl_tr=[]
    ensembl_pr=[]
    uniprot_id=[]
    uniprot_entry=[]
    refseq=[]
    for key in psm_hash:
        if 'search_hit' in key:
            for psm in key['search_hit']:
                for i in range(0, len(psm['proteins'])):
                    find_result=_find_annotation_individually_(psm['proteins'][i]['protein'],species)
                    if find_result!=None:
                        psm['proteins'][i]['protein'] = _update_protein_accession_(psm['proteins'][i]['protein'],
                                                                                   decoy_annotation, find_result[2])
                        hit = 0
                        if 'DECOY_' in psm['proteins'][i]['protein'].upper():
                            hit = 1
                        if hit == 0:
                            if find_result[1]=="ENSEMBL_TR":
                                ensembl_tr.append(psm['proteins'][i]['protein'])
                            elif find_result[1]=="ENSEMBL_PR":
                                ensembl_pr.append(psm['proteins'][i]['protein'])
                            elif find_result[1]=="UNIPROT_ENTRY":
                                uniprot_entry.append(psm['proteins'][i]['protein'])
                            elif find_result[1]=="UNIPROT":
                                uniprot_id.append(psm['proteins'][i]['protein'])
                            elif find_result[1]=="REFSEQ":
                                refseq.append(psm['proteins'][i]['protein'])
    if refseq==[] and uniprot_id==[] and uniprot_entry==[] and ensembl_pr==[] and ensembl_tr==[]:
        raise ValueError('Protein ID annotations not recognized\n' \
                         'currently supported Protein ID annotations:\n' \
                         'ENSEMBL protein identifiers\n ENSEMBL transcript identifiers\n UNIPROT/SWISSPROT identifiers '
                         '\n UNIPROT ENTRIES \n REFSEQ identifiers')
        sys.exit()
    id_map={}
    temp_psm_hash={}
    temp_exon_hash={}
    annotation=[temp_psm_hash,temp_exon_hash]
    if ensembl_tr!=[]:
        print "Identified ENSEMBL transcript IDs"
        id_map.update(_id_map_('ENSEMBL', 'ENSEMBL', ensembl_tr, psm_hash, species, decoy_annotation, database_v))
        temp_annotation = proBAM_ENSEMBL.prepareAnnotationENSEMBL(_get_transcript_ids_from_map(id_map), 'transcript',
                                                             database_v, species, three_frame_translation)
        temp_psm_hash.update(temp_annotation[0])
        temp_exon_hash.update(temp_annotation[1])
    if ensembl_pr!=[]:
        print "Identified ENSEMBL protein IDs"
        id_map.update(_id_map_('ENSEMBL', 'ENSEMBL', ensembl_pr, psm_hash, species, decoy_annotation, database_v))
        temp_annotation = proBAM_ENSEMBL.prepareAnnotationENSEMBL(_get_transcript_ids_from_map(id_map), 'protein',
                                                             database_v, species, three_frame_translation)
        temp_psm_hash.update(temp_annotation[0])
        temp_exon_hash.update(temp_annotation[1])
    if uniprot_id!=[]:
        temp_id_map_uniprot=_id_map_('UNIPROT', 'ENSEMBL', uniprot_id, psm_hash, species, decoy_annotation, database_v)
        if len(temp_id_map_uniprot.keys())>0:
            print " Identified UNIPROT IDs"
            id_map.update(temp_id_map_uniprot)
            temp_annotation = proBAM_ENSEMBL.prepareAnnotationENSEMBL(_get_transcript_ids_from_map(id_map), 'transcript',
                                                                 database_v, species, three_frame_translation)
            temp_psm_hash.update(temp_annotation[0])
            temp_exon_hash.update(temp_annotation[1])
    if uniprot_entry!=[]:
        print " Identified UNIPROT Entries"
        id_map.update(_id_map_('UNIPROT_ENTRY', 'ENSEMBL', uniprot_entry, psm_hash, species, decoy_annotation,
                          database_v))
        temp_annotation = proBAM_ENSEMBL.prepareAnnotationENSEMBL(_get_transcript_ids_from_map(id_map), 'transcript',
                                                             database_v, species, three_frame_translation, )
        temp_psm_hash.update(temp_annotation[0])
        temp_exon_hash.update(temp_annotation[1])
    if refseq!=[]:
        print " Identified RefSeq ID"
        id_map.update(_id_map_('REFSEQ', 'ENSEMBL', uniprot_entry, psm_hash, species, decoy_annotation,
                          database_v))
        temp_annotation = proBAM_ENSEMBL.prepareAnnotationENSEMBL(_get_transcript_ids_from_map(id_map), 'transcript',
                                                             database_v, species, three_frame_translation, )
        temp_psm_hash.update(temp_annotation[0])
        temp_exon_hash.update(temp_annotation[1])
    return [psm_hash, annotation, id_map]
#
# Create a list of unique transcript ID's from MAP
#
def _get_transcript_ids_from_map(map):
    '''
    :param map: dictionairy with protein IDs
    :return: unique list of transcript ids
    '''
    transcript_ids=[]
    for value in map.values():
        transcript_ids+=value
    transcript_ids=list(set(transcript_ids))
    return transcript_ids

#
# Find annotation without print (use: for single ID retrieval)
#
def _find_annotation_individually_(accession,species):
    '''
        :param term_list: list of terms of possible identifiers
        :param species: species
        :return: list [found_boolean,Type_identified_identifier,position of identifier term]
        '''
    ensembl_prefix = proBAM_ENSEMBL.get_Ensembl_prefix(species)
    result=None
    if re.findall(ensembl_prefix[0][0] + "[0-9]*", accession) != []:
        result = [1, "ENSEMBL_TR", ensembl_prefix[0][0] + "[0-9]*"]
    elif re.findall(ensembl_prefix[0][1] + "[0-9]*", accession) != []:
        result = [1, "ENSEMBL_PR", ensembl_prefix[0][1] + "[0-9]*"]
    elif re.findall("[A-Z0-9]{1,10}" + "_" + _get_uniprot_postfix_(species), accession) != []:
        result = [1, "UNIPROT_ENTRY", "[A-Z0-9]{1,10}" + "_" + _get_uniprot_postfix_(species)]
    elif re.findall("[OPQ][0-9][A-Z0-9]{3}[0-9]|[A-NR-Z][0-9][A-Z][A-Z0-9]{2}[0-9]|[A-NR-Z][0-9][A-Z][A-Z0-9]{2}[0-9][A-Z][A-Z0-9]{2}[0-9]", accession) != []:
        result = [1, "UNIPROT",
                  "[OPQ][0-9][A-Z0-9]{3}[0-9]|[A-NR-Z][0-9][A-Z][A-Z0-9]{2}[0-9]|[A-NR-Z][0-9][A-Z][A-Z0-9]{2}[0-9][A-Z][A-Z0-9]{2}[0-9]"]
    elif re.findall("[NX][MP][_][0-9]{9}|[NX][R][_][0-9]{6}", accession)!=[]:
        result = [1, "REFSEQ","[NX][MP][_][0-9]{9}|[NX][R][_][0-9]{6}"]
    return result
#
# Find annotation in list of terms
#
def _find_annotation_(accession,species,pre_picked_annotation):
    '''
    :param term_list: list of terms of possible identifiers
    :param species: species
    :return: list [found_boolean,Type_identified_identifier,position of identifier term]
    '''
    found=0
    ensembl_prefix=proBAM_ENSEMBL.get_Ensembl_prefix(species)
    result=[0,'']
    if re.findall(ensembl_prefix[0][0]+"[0-9]*",accession)!=[] and \
            (pre_picked_annotation=="First" or pre_picked_annotation=="Ensembl_tr") and found==0:
        print "Identified ENSEMBL transcript IDs"
        result=[1,"ENSEMBL_TR",ensembl_prefix[0][0]+"[0-9]*"]
        found=1
    if re.findall(ensembl_prefix[0][1]+"[0-9]*",accession)!=[] and \
        (pre_picked_annotation == "First" or pre_picked_annotation == "Ensembl_pr") and found==0 :
        print "Identified ENSEMBL protein IDs"
        result=[1,"ENSEMBL_PR",ensembl_prefix[0][1]+"[0-9]*"]
        found=1
    if (pre_picked_annotation=="First" or pre_picked_annotation=="UniProt_ACC") and found==0:
        if re.findall("[OPQ][0-9][A-Z0-9]{3}[0-9]|[A-NR-Z][0-9][A-Z][A-Z0-9]{2}[0-9]|[A-NR-Z][0-9][A-Z][A-Z0-9]{2}[0-9][A-Z][A-Z0-9]{2}[0-9]",accession)!=[]:
            print " Identified UNIPROT IDs"
            result=[1,"UNIPROT","[OPQ][0-9][A-Z0-9]{3}[0-9]|[A-NR-Z][0-9][A-Z][A-Z0-9]{2}[0-9]|[A-NR-Z][0-9][A-Z][A-Z0-9]{2}[0-9][A-Z][A-Z0-9]{2}[0-9]"]
            found=1
    if (pre_picked_annotation == "First" or pre_picked_annotation == "UniProt_Entry") and found == 0:
        if re.findall("[A-Z0-9]{1,10}"+"_"+_get_uniprot_postfix_(species),accession)!=[]:
            print "Identified UNIPROT Entry Names"
            result=[1,"UNIPROT_ENTRY","[A-Z0-9]{1,10}"+"_"+_get_uniprot_postfix_(species)]
            found=1
    if (pre_picked_annotation == "First" or pre_picked_annotation == "RefSeq") and found == 0:
        if re.findall("[NX][MP][_][0-9]{9}|[NX][R][_][0-9]{6}",accession)!=[]:
            print "Identified RefSeq IDs"
            result=[1,"REFSEQ","[NX][MP][_][0-9]{9}|[NX][R][_][0-9]{6}"]
            found=1
    return result



#
# Parse protein accesions
#

def _update_protein_accession_(accession,decoy_annotation,regex):
    '''
    :param accession: list of possible protein ID terms
    :param decoy_annotation: list of decoy annotations
    :param hit: position of protein ID
    :return: protein ID
    '''
    is_decoy=0

    if any(decoy in accession.upper().replace('-','_') for decoy in decoy_annotation):
        is_decoy=1
    new_accession=re.findall(regex,accession)
    if new_accession!=[] and new_accession[0]!="P00000":
        accession= new_accession[0]
        if is_decoy==1:
            accession="DECOY_"+accession
    return accession
#
# Get UNIPROT entry names postfixes
#
def _get_uniprot_postfix_(species):
    '''
    :param species: species
    :return: species dependent postfix for UniProt Entry ID's
    '''
    if species=='homo_sapiens':
        return "HUMAN"
    elif species=='mus_musculus':
        return "MOUSE"
    elif species=='drosophila_melanogaster':
        return "DROME"
    elif species=='danio_rerio':
        return "DANRE"
    elif species=='arabidopsis_thaliana':
        return "ARATH"
    else:
        raise ValueError('Species not recognized')
#
# Retrieve matched protein IDs
#
def _get_protein_ID_(psm_hash,decoy_annotation):
    '''
    :param psm_hash: dictionairy containing PSM file
    :param decoy_annotation: list of decoy annotation tags
    :return: list of unique protein IDs
    '''

    #retrieve protein IDs from PSM HASH
    protein_ID=[]
    for key in psm_hash:
        if 'search_hit' in key:
            for psm in key['search_hit']:
                for i in range (0,len(psm['proteins'])):
                    hit=0
                    for d in decoy_annotation:
                        if d in psm['proteins'][i]['protein'].upper():
                            ID=psm['proteins'][i]['protein'].upper().strip(d)
                            protein_ID.append(ID)
                            hit=1
                    if hit==0:
                        protein_ID.append(psm['proteins'][i]['protein'])
    #retain only unique IDs
    unique=set(protein_ID)
    #print unique
    return list(unique)

#
# Function to map protein IDs between annotations
#
def _id_map_(from_annotation,to_annotation,psm_protein_id,psm_hash,species,decoy_annotation,database_v):
    '''
    :param from_annotation: supplied annotation (i.e. swissprot)
    :param to_annotation: target annotation (i.e. ENSEMBL)
    :param psm_protein_id: list of protein IDS
    :param psm_hash: dictionairy of protein IDs mapped onto ENSEMBL
    :param species: species name
    :param decoy_annotation: list of decoy annotations
    :param database_v: database version
    :return: dictionairy of protein ID coversion
    '''
    #psm_hash.reset()
    new_psm_protein_id=[]
    psm_protein_id=list(set(psm_protein_id))
    print "Commencing ID conversion from "+str(from_annotation)+" to "+str(to_annotation)
    map={}

    # Convert RefSeq to ENSEMBL ID's
    if to_annotation=="ENSEMBL":
        if from_annotation=='REFSEQ':
            temp_map={}
            refseq_mrna=[]
            refseq_pred_mrna=[]
            refseq_ncrna=[]
            refseq_pred_ncrna=[]
            refseq_prot=[]
            refseq_pred_prot=[]

            for id in psm_protein_id:
                if "NM_" in id:
                    refseq_mrna.append(id)
                if "XM_" in id:
                    refseq_pred_mrna.append(id)
                if "NR_" in id:
                    refseq_ncrna.append(id)
                if "XR_" in id:
                    refseq_pred_ncrna.append(id)
                if "NP_" in id:
                    refseq_prot.append(id)
                if "XP_" in id:
                    refseq_pred_prot.append(id)
            mapped_id=[]
            if refseq_mrna!=[]:
                print "Identified refseq mRNA ID's, converting:"
                mapped_id.append(proBAM_biomart.id_map_ensembl("refseq_mrna", database_v, species, refseq_mrna))
            if refseq_pred_mrna!=[]:
                print "Identified predicted refseq mRNA ID's, converting:"
                mapped_id.append(proBAM_biomart.id_map_ensembl("refseq_mrna_predicted", database_v, species, refseq_pred_mrna))
            if refseq_ncrna!=[]:
                print "Identified refseq ncRNA ID's, converting:"
                mapped_id.append(proBAM_biomart.id_map_ensembl("refseq_ncrna", database_v, species, refseq_ncrna))
            if refseq_pred_ncrna!=[]:
                print "Identified refseq predicted ncRNA ID's, converting:"
                mapped_id.append(proBAM_biomart.id_map_ensembl("refseq_ncrna_predicted", database_v, species, refseq_pred_ncrna))
            if refseq_prot!=[]:
                print "Identified refseq protein ID's, converting:"
                mapped_id.append(proBAM_biomart.id_map_ensembl("refseq_peptide", database_v, species, refseq_prot))
            if refseq_pred_prot!=[]:
                print "Identified refseq predicted protein ID's, converting:"
                mapped_id.append(proBAM_biomart.id_map_ensembl("refseq_peptide_predicted", database_v, species, refseq_pred_prot))
            for row in mapped_id:
                if row[0]!="":
                    if row[2] in temp_map:
                            temp_map[row[2]].append(row)
                    else:
                        temp_map[row[2]]=[row]
            for key in temp_map:
                map[key]=temp_map[key][0]

        #Convert UNIPROT accession ID's to ENSEMBL
        if from_annotation=="UNIPROT":
            temp_map={}
            #map uniprot/swissprot
            mapped_id=proBAM_biomart.id_map_ensembl("uniprot_swissprot",database_v,species,psm_protein_id)
            if len(mapped_id)>1:
                for row in mapped_id:
                    if row[0]!="":
                        if row[2] in temp_map:
                                temp_map[row[2]].append(row)
                        else:
                            temp_map[row[2]]=[row]
                for key in temp_map:
                    map[key]=temp_map[key][0]

            #map remaining on uniprot/trembl
            unmapped_id_for_trmbl=[]
            for id in psm_protein_id:
                if id not in map:
                    unmapped_id_for_trmbl.append(id)
            if unmapped_id_for_trmbl!=[]:
                mapped_id=proBAM_biomart.id_map_ensembl("uniprot_sptrembl",database_v,species,unmapped_id_for_trmbl)
                if len(mapped_id)>1:
                    for row in mapped_id:
                        if row[0]!="":
                            if row[2] in temp_map:
                                    temp_map[row[2]].append(row)
                            else:
                                temp_map[row[2]]=[row]
                    for key in temp_map:
                        map[key]=temp_map[key][0]

        #Convert UNIPROT ENTRIES to ENSEMBL ID's
        if from_annotation=="UNIPROT_ENTRY":
            u=UniProt()
            to_translate=[]
            for id in psm_protein_id:
                if re.findall("[A-Z0-9]{1,10}" + "_" + _get_uniprot_postfix_(species), id) != []:
                    to_translate.append(re.findall("[A-Z0-9]{1,10}" + "_" + _get_uniprot_postfix_(species), id)[0])
            to_translate = list(set(to_translate))
            if len(to_translate) > 1000:
                nr_chunks = len(to_translate) / 1000
                chunks = chunkIt(to_translate, nr_chunks)
            else:
                chunks = [to_translate]
            # map uniprot_entries to up-to-date accession
            accession_update_hash={}
            for chunk in chunks:
                chunk_accession_update_hash=u.mapping('ACC+ID','ACC',chunk)
                accession_update_hash.update(chunk_accession_update_hash)

            tot_count=0
            for entry in to_translate:
                if entry not in accession_update_hash:
                    tot_count+=1

            print "\tFound "+str(tot_count)+" depreciated UniProt Entries (UniProt Entries are unstable)\n" \
                  "\tAttempting to map depreciated entries unto new entry ID's..."
            count=0
            found=0
            tracker=9.99
            print "\t",
            for entry in to_translate:
                if entry not in accession_update_hash:
                    count+=1
                    if (float(count * 100) / float(tot_count))>tracker:
                        print str(int(tracker+0.01))+"% ",
                        tracker=tracker+10
                    try:
                        accession_update_hash[entry]=[get_updated_entry_name(entry)]
                        found+=1
                    except urllib2.HTTPError:
                        pass
            print ' '
            to_translate=[]
            print "\tRetrieved "+str(found)+" of the "+str(tot_count)+" depreciated UniProt Entries"

            #map accession to Ensembl transcript
            for value in accession_update_hash.values():
                to_translate+=value
            to_translate = list(set(to_translate))
            if len(to_translate) > 1000:
                nr_chunks = len(to_translate) / 1000
                chunks = chunkIt(to_translate, nr_chunks)
            else:
                chunks = [to_translate]
            # map uniprot_entries to up-to-date accession
            temp_map={}
            for chunk in chunks:
                # remap Ensembl transcript to uniprot entries
                chunk_temp_map=u.mapping('ACC','ENSEMBL_TRS_ID',chunk)
                temp_map.update(chunk_temp_map)
            for accession in accession_update_hash:
                for i in accession_update_hash[accession]:
                    if i in temp_map:
                        if accession in map:
                            map[accession]+=temp_map[i]
                        else:
                            map[accession]=temp_map[i]
        #Get ENSEMBL ID's in correct form
        if from_annotation=="ENSEMBL":
            for id in psm_protein_id:
                map[id]=[id]
    return map
#
# Divides a array in equal parts (chunks) ( fair redistributing of computational load)
#
def chunkIt(seq, num):
    '''
    :param seq: sequence or array
    :param num: number of chunks to be distributed
    :return: an array with in each compartment an equal amount of values
    '''
    avg = len(seq) / float(num)
    out = []
    last = 0.0

    while last < len(seq):
        out.append(seq[int(last):int(last + avg)])
        last += avg

    return out
#
# given a UniProt entry name, retrieves the stable accession ID
#
def get_updated_entry_name(name):
    '''
    :param name: UniProt Entry
    :return: uniprot accession
    '''
    url = 'http://www.uniprot.org/uniprot/'+name+'.tab'

    params = {
    'columns':'id',
    }

    data = urllib.urlencode(params)
    request = urllib2.Request(url, data)
    contact = "volodimir.olexiouk@ugent.be" # Please set your email address here to help us debug in case of problems.
    request.add_header('User-Agent', 'Python %s' % contact)
    response = urllib2.urlopen(request)
    page = response.read(200000)
    return page.split("\n")[1]

