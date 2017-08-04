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

from cogent.core.genetic_code import DEFAULT as standard_code
import proBAM_biomart
import MySQLdb
import time


#
# Retrieve correct Ensembl identifiers
#
def get_Ensembl_prefix(species):
    '''
    :param species: fulle species name
    :return: list of ENSEMBL transcript and protein tags
    '''
    ensembl_prefix=[]
    if species=='homo_sapiens':
        ensembl_prefix.append(['ENST','ENSP'])
    elif species=='mus_musculus':
        ensembl_prefix.append(['ENSMUST','ENSMUSP'])
    elif species=='drosophila_melanogaster':
        ensembl_prefix.append(['FBtr','FBpp'])
    elif species=='danio_rerio':
        ensembl_prefix.append(['ENSDART','ENSDARP'])
    elif species=='arabidopsis_thaliana':
        ensembl_prefix.append(['AT','AT'])
    else:
        raise ValueError('Species not recognized')
    return ensembl_prefix

#
# Retrieve data from BioMart and store in hashes
#
def prepareAnnotationENSEMBL(psm_protein_id,mode,database_v,species,three_frame_translation):
    '''
    :param psm_protein_id: list of protein IDs (untagged)
    :param mode: transcript or protein mode
    :param database_v: database version
    :param species: species name
    :return: dictionairy mapping proteins into ENSEMBL
    '''
    try:
        db = MySQLdb.connect(host='ensembldb.ensembl.org',user='anonymous',passwd='',port=3306)
        cur = db.cursor()
    except:
        print "Could not connect to Ensembl trough port 3306, attempting to connect trough port 5306"
        db = MySQLdb.connect(host='ensembldb.ensembl.org', user='anonymous', passwd='', port=5306)
        cur = db.cursor()
        pass

    cur.execute('show databases')
    for row in cur.fetchall():
        if species+'_core_'+str(database_v) in row[0]:
            mysql_database=row[0]
            break
    db.close()
    try:
        db = MySQLdb.connect(host='ensembldb.ensembl.org',user='anonymous',passwd='',port=3306,db=mysql_database)
        cur = db.cursor()
    except:
        try:
            print "Could not connect to Ensembl trough port 3306, attempting to connect trough port 5306"
            db = MySQLdb.connect(host='ensembldb.ensembl.org', user='anonymous', passwd='', port=5306,db=mysql_database)
            cur = db.cursor()
        except:
            print "could not connect to Ensembl, check firewall/internet connection"
            raise


    chunked_psm_protein_id=chunkIt(psm_protein_id,10)
    process=0
    psm_protein_id = {}
    transcript_ids = []
    for chunk in chunked_psm_protein_id:
        first=1
        where_clause='WHERE '
        for i in range(0,len(chunk)):
            if mode == 'protein':
                if first==1:
                    where_clause+="translation.stable_id='"+str(chunk[i])+"' "
                    first=0
                else:
                    where_clause += "OR translation.stable_id='" + str(chunk[i]) + "' "
            elif mode == 'transcript':
                if first==1:
                    where_clause+="transcript.stable_id='"+str(chunk[i])+"' "
                    first=0
                else:
                    where_clause += "OR transcript.stable_id='" + str(chunk[i]) + "' "

        if mode == 'protein':
            id = 1
            sql='SELECT transcript.stable_id,translation.stable_id,transcript.transcript_id,'\
                        'translation.seq_start,translation.start_exon_id '\
                        'from transcript LEFT JOIN translation ON '\
                        'transcript.transcript_id=translation.transcript_id '+where_clause
            cur.execute(sql)

        elif mode == 'transcript':
            id = 0
            sql='SELECT transcript.stable_id,translation.stable_id,transcript.transcript_id,'\
                        'translation.seq_start,translation.start_exon_id '\
                        'FROM transcript LEFT JOIN translation ON '\
                        'transcript.transcript_id=translation.transcript_id '+where_clause
            cur.execute(sql)

        for row in cur.fetchall():

            transcript_ids.append(row[2])
            if row[1]!=None or three_frame_translation=='Y':
                psm_protein_id[row[id]] = {'transcript_id': row[0], 'translation_id': row[1],
                                           'transcript_seq': '', 'protein_seq': '',
                                           'chr': '', 'strand': '', '5UTR_offset': row[3], 'start_exon_rank': row[4]}

        if process < 100:
            process += 10
            print str(process) + "% ",

    print " "
    return ensembl_construct_sequences(psm_protein_id,mysql_database,transcript_ids,database_v,species,
                                       three_frame_translation,mode)

#
# Retrieve DNA-sequence, Splice, construct protein sequence and store
#
def ensembl_construct_sequences(psm_hash,mysql_db,transcript_ids,database_v,species,three_frame_translation,mode,):
    '''
    :param psm_hash: dictionair with protein / ensembl information ( see prepareAnnotationENSEMBL)
    :param ensembl:ensembl genome
    :param transcript_ids: list of transcrip ids (converted from protein IDs)
    :param database_v: database version
    :param species: species name
    :return: dictionairy mapping proteins into ENSEMBL
    '''
    print "Commencing transcript and protein sequence retrieval"

    no_protein_seq=[]
    biomart_key_hash={}
    stable_transcript_ids=[]

    for key in psm_hash.keys():
        biomart_key_hash[psm_hash[key]['transcript_id']]=key
        stable_transcript_ids.append(psm_hash[key]['transcript_id'])

    chunked_stable_transcript_id=chunkIt(stable_transcript_ids,10)
    process=0
    c=0
    for chunk in chunked_stable_transcript_id:
        # Retrieve cds,chr,transcript_id and strand from biomart
        try:
            biomart_result=proBAM_biomart.retrieve_data_from_biomart(database_v,species,chunk,three_frame_translation)
        except AttributeError:
            time.sleep(60)
            print "BioMart connection timeout, reconnecting to BioMart"
            biomart_result = proBAM_biomart.retrieve_data_from_biomart(database_v, species, chunk, three_frame_translation)
        for row in biomart_result:
            row=row.split("\t")
            try:
                psm_hash[biomart_key_hash[row[1]]]['transcript_seq']=row[0]
                psm_hash[biomart_key_hash[row[1]]]['shift']=_calc_seq_shift_(row[0])
                psm_hash[biomart_key_hash[row[1]]]['protein_seq']=standard_code.translate(row[0])
                psm_hash[biomart_key_hash[row[1]]]['chr']=row[2]
                psm_hash[biomart_key_hash[row[1]]]['strand']=row[3]
                del row
            except IndexError:
                pass
        del biomart_result
        if process<100:
            process+=10
            print str(process)+"% ",
    print " "

    # get exons directly from core database
    temp_exon_hash=get_ensembl_exons(mysql_db,transcript_ids,psm_hash,mode)
    exon_hash=temp_exon_hash[0]
    psm_hash=temp_exon_hash[1]
    del temp_exon_hash
    # retrieve protein sequences for transcript where the protein sequence could not be fetched automatically
    for key in no_protein_seq:
        psm_hash[key]['transcript_seq']=retrieve_protein_seq(psm_hash[key]['transcript_seq'],
                                                                   exon_hash[psm_hash[key]['transcript_id']],
                                                                   psm_hash[key]['5UTR_offset'],
                                                                   psm_hash[key]['start_exon_rank'])
        psm_hash[key]['shift'] = _calc_seq_shift_(psm_hash[key]['transcript_seq'])
        #translate till stop codon
        psm_hash[key]['protein_seq']=standard_code.translate(psm_hash[key]
                                                                   ['transcript_seq']).partition('*')[0]
    return [psm_hash,exon_hash]

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
#
#
def _calc_seq_shift_(sequence):

    i = 0
    hit=0
    shift = 0
    while hit==0:
        if sequence[i] == 'N':
            shift += 1
        else:
            hit = 1
        i += 1
    return shift
#
# Retrieve protein coding sequence from transcript
#

def retrieve_protein_seq(transcript_seq,exons,offset,start_exon):
    '''
    :param transcript_seq: transcript sequence
    :param exons: exons from transcript sequence
    :param offset: offset from transcript start site till CDS start site
    :param start_exon: first exon containing CDS
    :return: CDS sequence
    '''
    count=0
    for exon in exons:
        exon=[int(numeric_string) for numeric_string in exon]
        if exon[2]<start_exon:
            count+=exon[1]-exon[0]
        elif exon[2]==start_exon:
            count+=offset
            break
    return transcript_seq[(count-1):len(transcript_seq)]

def get_ensembl_exons(mysql_database,transcript_ids,psm_hash,mode):
    '''
    :param ensembl: ENSEMBL genome
    :param transcript_ids: list of transcript ids
    :param psm_hash: dictionairy mapping proteins as transcript on ensembl
    :return: exon hash
    '''
    print "Commencing exon retrieval"
    transcript_ids=list(set(transcript_ids))
    exon_hash={}
    prot_tr = {}
    #connect and join tables in order to retrieve exon information
    chunked_transcript_ids=chunkIt(transcript_ids, 10)
    process=0
    try:
        db = MySQLdb.connect(host='ensembldb.ensembl.org',user='anonymous',passwd='',port=3306, db=mysql_database)
        cur = db.cursor()
    except:
        try:
            print "Could not connect to Ensembl trough port 3306, attempting to connect trough port 5306"
            db = MySQLdb.connect(host='ensembldb.ensembl.org', user='anonymous', passwd='', port=5306, db=mysql_database)
            cur = db.cursor()
        except:
            print "could not connect to Ensembl, check firewall/internet connection"
            raise
    for chunk in chunked_transcript_ids:
        first=1
        where_clause='WHERE '
        for i in range(0,len(chunk)):
            if first == 1:
                where_clause += "exon_transcript.transcript_id='" + str(chunk[i]) + "' "
                first = 0
            else:
                where_clause += "OR exon_transcript.transcript_id='" + str(chunk[i]) + "' "

        cur.execute("SELECT exon.seq_region_start, exon.seq_region_end, exon_transcript.rank, transcript.stable_id, " \
        "exon_transcript.exon_id " \
        "FROM exon_transcript LEFT JOIN exon " \
        "ON exon.exon_id=exon_transcript.exon_id " \
        "LEFT JOIN transcript " \
        "ON transcript.transcript_id=exon_transcript.transcript_id " + where_clause)

        #speed up process by making a hash
        if mode=='protein':
            for key in psm_hash:
                prot_tr[psm_hash[key]['transcript_id']]=key
        for row in cur.fetchall():
            if row[3] not in exon_hash:
                exon_hash[row[3]]=[]
            if mode=='transcript':
                if row[3] in psm_hash:
                    if psm_hash[row[3]]['start_exon_rank']==row[4]:
                        psm_hash[row[3]]['start_exon_rank']=row[2]
            elif mode=='protein':
                if row[3] in prot_tr:
                    if psm_hash[prot_tr[row[3]]]['start_exon_rank'] == row[4]:
                        psm_hash[prot_tr[row[3]]]['start_exon_rank'] = row[2]
            exon_hash[row[3]].append([str(row[0]),str(row[1]),str(row[2])])
        if process < 100:
            process += 10
            print str(process) + "% ",

    db.close()
    print " "
    return [exon_hash,psm_hash]

#
# Create SQ header for SAM file (ENSEMBL) get chr location and coordinates and store them in an array
#
def create_SQ_header(database_v,species):
    '''
    :param database_v: database version
    :return: list of chromosomes with their size (from ENSEMBL)
    '''
    SQ=[]
    # create connection to ensembm database
    # create connection to ensembl database
    try:
        db = MySQLdb.connect(host='ensembldb.ensembl.org',user='anonymous',passwd='',port=3306)
        cur = db.cursor()
    except:
        try:
            print "Could not connect to Ensembl trough port 3306, attempting to connect trough port 5306"
            db = MySQLdb.connect(host='ensembldb.ensembl.org', user='anonymous', passwd='', port=5306)
            cur = db.cursor()
        except:
            print "could not connect to Ensembl, check firewall/internet connection"
            raise

    cur.execute('show databases')
    for row in cur.fetchall():
        if species + '_core_' + str(database_v) in row[0]:
            mysql_database = row[0]

    db.close()
    try:
        db = MySQLdb.connect(host='ensembldb.ensembl.org',user='anonymous',passwd='',port=3306, db=mysql_database)
        cur = db.cursor()
    except:
        try:
            print "Could not connect to Ensembl trough port 3306, attempting to connect trough port 5306"
            db = MySQLdb.connect(host='ensembldb.ensembl.org', user='anonymous', passwd='', port=5306, db=mysql_database)
            cur = db.cursor()
        except:
            print "could not connect to Ensembl, check firewall/internet connection"
            raise
    cur.execute("select seq_region.name,seq_region.length,coord_system.version FROM seq_region LEFT JOIN coord_system "
                "ON coord_system.coord_system_id=seq_region.coord_system_id WHERE coord_system.rank=1")
    for row in cur.fetchall():
        if '_' not in row[0]:
            SQ_string= "@SQ\tSN:chr"+str(row[0])+"\tLN:"+str(row[1])+"\tAS:"+str(row[2])+"\tSP:"+str(species)
            SQ.append(SQ_string)
    db.close()
    return SQ

#
# get Genome assembly ID
#
#create_SQ_headeron and coordinates and store them in an array
#
def get_genome_version(database_v,species):
    '''
    :param database_v: database version
    :return: list of chromosomes with their size (from ENSEMBL)
    '''
    SQ=[]
    # create connection to ensembm database
        # create connection to ensembl database
    try:
        db = MySQLdb.connect(host='ensembldb.ensembl.org',user='anonymous',passwd='',port=3306)
        cur = db.cursor()
    except:
        try:
            print "Could not connect to Ensembl trough port 3306, attempting to connect trough port 5306"
            db = MySQLdb.connect(host='ensembldb.ensembl.org', user='anonymous', passwd='', port=5306)
            cur = db.cursor()
        except:
            print "could not connect to Ensembl, check firewall/internet connection"
            raise

    cur.execute('show databases')
    for row in cur.fetchall():
        if species+'_core_'+str(database_v) in row[0]:
            mysql_database=row[0]
            break
    db.close()
    try:
        db = MySQLdb.connect(host='ensembldb.ensembl.org',user='anonymous',passwd='',port=3306, db=mysql_database)
        cur = db.cursor()
    except:
        try:
            print "Could not connect to Ensembl trough port 3306, attempting to connect trough port 5306"
            db = MySQLdb.connect(host='ensembldb.ensembl.org', user='anonymous', passwd='', port=5306, db=mysql_database)
            cur = db.cursor()
        except:
            print "could not connect to Ensembl, check firewall/internet connection"
            raise

    cur.execute("SELECT version from coord_system where rank=1")

    for row in cur.fetchall():
        result=row[0]

    version=str(species+"."+result)
    db.close()
    return version
