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
from pyteomics import pepxml
from fnmatch import fnmatch
from lxml import etree


#
# Reads pepXML and stores variables in an hash
#
def get_PSM_pepxml(psm_file):
    '''
    :param psm_file: psm file
    :return: dictionairy of psms
    '''
    PSM = []
    PEP = pepxml.read(psm_file, read_schema=False, iterative=True)
    # count = 0
    # parse tags out of protein IDs
    for row in PEP:
        # adjust search scores
        if 'search_hit' in row.keys():
            for search_hit in row['search_hit']:
                search_hit['massdiff'] = str(
                    search_hit['massdiff']) + ';' + str(
                        search_hit['calc_neutral_pep_mass']) + ';' + str(
                            row['precursor_neutral_mass'])
                search_hit['search_score'] = {
                    'score': _get_score_(search_hit['search_score']),
                    'evalue': _get_evalue_(search_hit['search_score'])
                }

        # if count==5500:
        #    break
        PSM.append(row)
        # count+=1
    del PEP
    return PSM


#
# Attempts to extract enzyme information from the pepXML file
#
def get_enzyme_pepxml(psm_file):
    '''
    :param psm_file: psm_file
    :return: enzyme
    '''
    try:
        root = etree.parse(psm_file).getroot()
        namespace = root.nsmap.items()[0][1]
        line = root.xpath(
            "//ns:sample_enzyme", namespaces={'ns': namespace})[0].get("name")
        line = line.replace('-', '_')
        line = line.lower()
        if fnmatch(line, '*trypsin?p*'):
            return 2
        elif fnmatch(line, '*trypsin*'):
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
        elif fnmatch(line, '*no?enzyme*'):
            return 0
        return "-1"
    except:
        return "-1"


# Attempts to extract enzyme specificity information from the pepXML file
#
def get_enzyme_specificity_pepxml(psm_file):
    '''
    :param psm_file: psm file
    :return: enzyme specificity
    '''
    try:
        root = etree.parse(psm_file).getroot()
        namespace = root.nsmap.items()[0][1]
        line = root.xpath(
            "//ns:specificity", namespaces={'ns': namespace})[0].get("cut")
        line = line.replace('-', '_')
        line = line.lower()
        if fnmatch(line, '*non*enzy*'):
            return 0
        elif fnmatch(line, '*semi*enzy*'):
            return 1
        elif fnmatch(line, '*fully*enzy*'):
            return 2
        return 3
    except:
        return 3


#
# extract search summary information from pepXML
#
def extract_comments_from_pepxml(psm_file):
    '''
    :param psm_file: psm_file
    :return: comments array
    '''
    f = open(psm_file, 'r')
    count = 0
    comments = []
    while count < 200:
        count += 1
        line = f.readline()
        if "search_summary" in line:
            while "/search_summary" not in line:
                comments.append(line)
                line = f.readline()
            return comments
    return comments


#
# retrieve score for a psm
#
def _get_score_(search_score):
    '''
    :param search_score: score dictionairy
    :return: score
    '''
    hit = 0
    hit_key = ''
    score = {}
    # print psm.keys()
    for key in search_score:
        if "score" in key.lower():
            hit = 1
            hit_key = key
    if hit == 1:
        return search_score[hit_key]
    else:
        return "-1"


#
# retrieve the evalue for a psm
#
def _get_evalue_(search_score):
    '''
    :param search_score: score dictionairy
    :return: evalye
    '''

    hit = 0
    hit_key = ''
    score = {}
    # print psm.keys()
    for key in search_score.keys():
        if "xcorr" in key.lower() or 'expectation' in key.lower() or 'confidence' in key.lower() \
                or "e_value" in key.lower().replace("-", "_") or 'evalue' in key.lower():
            hit = 1
            hit_key = key
    if hit == 1:
        return search_score[hit_key]
    else:
        return "-1"
