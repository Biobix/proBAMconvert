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

from bioservices import BioMart

#
#Function that links the correct database archive with version number
#
def _get_ensembl_archive_(version,species):
    '''
    :param version: Ensembl version
    :return: ENSEMBL repository for a specific version
    '''
    version=int(version)
    if species=='arabidopsis_thaliana':
        return 'plants.ensembl.org'
    else:
        d={}
        d[89]="www.ensembl.org"
        d[88]="mar2017.archive.ensembl.org"
        d[87]="dec2016.archive.ensembl.org"
        d[86]="oct2016.archive.ensembl.org"
        d[85]="jul2016.archive.ensembl.org"
        d[84]="mar2016.archive.ensembl.org"
        d[83]="dec2015.archive.ensembl.org"
        d[82]="sep2015.archive.ensembl.org"
        d[81]="jul2015.archive.ensembl.org"
        d[80]="may2015.archive.ensembl.org"
        d[79]="mar2015.archive.ensembl.org"
        d[78]="dec2014.archive.ensembl.org"
        d[77]="oct2014.archive.ensembl.org"
        d[76]="aug2014.archive.ensembl.org"
        d[75]="feb2014.archive.ensembl.org"
        d[74]="dec2013.archive.ensembl.org"
        d[73]="sep2013.archive.ensembl.org"
        d[72]="jun2013.archive.ensembl.org"
        d[71]="apr2013.archive.ensembl.org"
        d[70]="jan2013.archive.ensembl.org"
        d[69]="oct2012.archive.ensembl.org"
        d[68]="jul2012.archive.ensembl.org"
        d[67]="may2012.archive.ensembl.org"
        d[66]="feb2012.archive.ensembl.org"
        d[65]="dec2011.archive.ensembl.org"
        d[64]="sep2011.archive.ensembl.org"
        d[63]="jun2011.archive.ensembl.org"
        d[62]="apr2011.archive.ensembl.org"
        d[61]="feb2011.archive.ensembl.org"
        d[60]="nov2010.archive.ensembl.org"
        d[59]="aug2010.archive.ensembl.org"
        d[58]="may2010.archive.ensembl.org"
        d[57]="mar2010.archive.ensembl.org"
        d[56]="sep2009.archive.ensembl.org"
        d[55]="jul2009.archive.ensembl.org"
        d[54]="may2009.archive.ensembl.org"
        if version in d:
            return d[version]
        else:
            raise ValueError('unsupported ensembl version')


#
# Function to map spesies to ENSEMBL dataset names
#
def _get_ensembl_dataset_(species):
    '''
    :param species: full species name
    :return: ensembl species name
    '''
    d={}
    d['homo_sapiens']='hsapiens_gene_ensembl'
    d['mus_musculus']='mmusculus_gene_ensembl'
    d['drosophila_melanogaster']='dmelanogaster_gene_ensembl'
    d['danio_rerio']='drerio_gene_ensembl'
    d['arabidopsis_thaliana']='athaliana_eg_gene'

    if species not in d:
        print 'Error: unsupported species'
        print 'Currently supported species:'
        print d.keys()
        raise ValueError('unsupported species')
    return d[species]
#
# Function to create XML readable transcript_id query string
#
def _id_in_xml_query_(transcipt_id):
    '''
    :param transcipt_id: list of transcrip IDs
    :return: XML readable transcript ID query string
    '''
    query=""
    for tr in transcipt_id:
        query+=(str(tr)+",")
    query=query[:-1]
    return query

#
# Function that retrieves cds,strand,chr and ensembl_transcript_id from BioMart
#
def retrieve_data_from_biomart(version,species,transcript_id,three_frame_translation):
    '''
    :param version: Database version
    :param species: Full species name
    :param transcript_id: list of transcript IDs
    :return: BioMart results
    '''

    #create connection
    tr_query=_id_in_xml_query_(transcript_id)
    version=_get_ensembl_archive_(version,species)

    dataset=_get_ensembl_dataset_(species)
    biomart = BioMart(host=version)

    #add filters
    biomart.add_dataset_to_xml(dataset)
    biomart.add_filter_to_xml("ensembl_transcript_id",tr_query)

    #add attributes
    biomart.add_attribute_to_xml('ensembl_transcript_id')
    biomart.add_attribute_to_xml("chromosome_name")
    biomart.add_attribute_to_xml("strand")
    if three_frame_translation=="Y":
        biomart.add_attribute_to_xml("cdna")
    else:
        biomart.add_attribute_to_xml("coding")
    attributes=biomart.attributes(dataset)

    #execute query
    xml_query=biomart.get_xml()

    # create bypass for plants database
    if species=="arabidopsis_thaliana":
        xml_query=xml_query.replace('virtualSchemaName = "default"','virtualSchemaName = "plants_mart_30"')

    result=biomart.query(xml_query)
    result=result.split("\n")

    return result

#
# Function that maps Identifiers to ENSEMBl
#
def id_map_ensembl(to_annotation,version,species,psm_protein_id):
    '''
    :param to_annotation: target identifier annotation (i.e. uniprot_swissprot)
    :param version: Database version
    :param species: Full species name
    :param psm_protein_id: list of IDs to be converted
    :return: BioMart results
    '''
    # If species is in plantsDB, execute plants adjusted function
    if species=="arabidopsis_thaliana":
        result=id_map_ensembl_plants(to_annotation,version,species,psm_protein_id)
        return result
    else:
        #adjust UniProt xml annotation for BioMart version >87
        if int(version)>87 and "uniprot" in to_annotation:
            to_annotation=to_annotation.replace('_','')
        #create connection
        query_string=_id_in_xml_query_(psm_protein_id)
        version=_get_ensembl_archive_(version,species)
        dataset=_get_ensembl_dataset_(species)
        biomart = BioMart(host=version)

        #add filters
        biomart.add_dataset_to_xml(dataset)
        biomart.add_filter_to_xml(to_annotation,query_string)

        #add attributs
        biomart.add_attribute_to_xml("ensembl_transcript_id")
        biomart.add_attribute_to_xml("transcript_start")
        biomart.add_attribute_to_xml("transcript_end")
        biomart.add_attribute_to_xml(to_annotation)
        attributes=biomart.attributes(dataset)

        #execute query
        xml_query=biomart.get_xml()
        tmp_result=biomart.query(xml_query)
        if len(tmp_result)==1:
            print "ERROR: could not convert ID's trough BioMart, " \
                  "Please check whether Ensembl version/species were correctly supplied"
        tmp_result=tmp_result.split("\n")
        result=[]

        if tmp_result!=[]:
            for i in tmp_result:
                i=i.split("\t")
                if i[0]!="":
                    result.append([i[0],(int(i[2])-int(i[1])),i[3]])
                else:
                    result.append(i)
        return result

#
# Function that maps Identfiers to Ensembl adjusted for plant DB compatibility
#

def id_map_ensembl_plants(to_annotation,version,species,psm_protein_id):
    '''
    :param to_annotation: to which annotation
    :param version: ensembl version
    :param species: species
    :param psm_protein_id: list of protein IDs
    :return: list of protein ID's converted to ENSEMBL
    '''
    #create connection
    query_string=_id_in_xml_query_(psm_protein_id)
    version=_get_ensembl_archive_(version,species)
    dataset=_get_ensembl_dataset_(species)
    biomart = BioMart(host=version)

    #add filters
    biomart.add_dataset_to_xml(dataset)
    biomart.add_filter_to_xml(to_annotation+"_accession",query_string)

    #add attributs
    biomart.add_attribute_to_xml("ensembl_transcript_id")
    biomart.add_attribute_to_xml("transcript_start")
    biomart.add_attribute_to_xml("uniprot_swissprot_accession")
    biomart.add_attribute_to_xml("transcript_end")

    #execute query
    xml_query=biomart.get_xml()
    xml_query=xml_query.replace('virtualSchemaName = "default"','virtualSchemaName = "plants_mart_30"')

    #parse results and adjust length
    temp_result=biomart.query(xml_query).split("\n")
    result=[]
    for row in temp_result:
        items=row.split("\t")
        # print row
        if len(items)==4:
            length=int(items[3])-int(items[1])+1
            result.append(items[0]+"\t"+str(length)+"\t"+items[2])
    return result

