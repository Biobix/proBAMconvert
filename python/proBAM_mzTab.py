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

import csv
from pyteomics import mass

#
# Parse mztab file into dictionairy
#


def get_PSM_mztab(psm_file):
    '''
    :param psm_file: mztab file
    :return: dictionairy of converted mztab file, suitable for proBAM conversion
    '''
    # connect to modification datbases and parse in desired format (UNIMOD & PSIMOD)
    unimod = _unimod_parser_()
    psimod = _psimod_xml_parser_()
    with open(psm_file, 'rU') as csvfile:
        #read mztab PSM rows
        mztab = csv.reader(csvfile, delimiter='\t')
        #hash to store variables in processable format
        psm_hash = []

        #iterate over all iterations and bundle psm from same spectrum
        spectrum = {}
        column_id = {}
        for row in mztab:
            psh_passed = False
            if row != "" and row != "\n" and len(row) != 0:
                if psh_passed == False and row[0] == "PSH":
                    for pos in range(0, len(row)):
                        column_id[row[pos]] = pos
                    for key in column_id.keys():
                        if "rank" in key.lower():
                            column_id['rank'] = int(column_id[key])
                        if "xcorr" in key.lower() or 'expectation' in key.lower() or 'confidence' in key.lower() \
                                or "e_value" in key.lower().replace("-","_") or 'evalue' in key.lower() or 'fdr' in key.lower():
                            column_id['fdr'] = int(column_id[key])
                    psh_passed = True

                if row[0] == "PSM":
                    if row[column_id['spectra_ref']] not in spectrum:
                        spectrum[row[column_id['spectra_ref']]] = [row]
                    else:
                        spectrum[row[column_id['spectra_ref']]].append(row)
        #iterate over all spectrum to store in processable dictionairy
        for key in spectrum.keys():
            if 'charge' in column_id:
                temp_hash = {
                    "assumed_charge": spectrum[key][0][column_id['charge']],
                    "spectrum": key,
                    "search_hit": []
                }
            else:
                temp_hash = {
                    "assumed_charge": 0,
                    "spectrum": key,
                    "search_hit": []
                }

            for psm in spectrum[key]:
                proteins = []
                proteins.append({
                    "protein": psm[column_id["accession"]],
                    'peptide_prev_aa': psm[column_id["pre"]],
                    "peptide_next_aa": psm[column_id["post"]]
                })
                modifications = _get_modifications_(
                    psm[column_id["modifications"]])
                #modified_sequence=_get_modified_sequence_(psm[column_id["sequence"]],psm[column_id["modifications"]],
                #                                          unimod,psimod)
                temp_hash['search_hit'].append({
                    "hit_rank":
                    _get_hit_rank_(psm, column_id),
                    "modifications":
                    modifications,
                    "calc_neutral_pep_mass":
                    psm[column_id['exp_mass_to_charge']],
                    "precursor_neutral_mass":
                    psm[column_id['calc_mass_to_charge']],
                    # "modified_peptide":modified_sequence,
                    "peptide":
                    psm[column_id['sequence']],
                    "massdiff":
                    _calc_massdiff_(psm[column_id['exp_mass_to_charge']],
                                    psm[column_id['calc_mass_to_charge']]),
                    "search_score": {
                        "score": psm[column_id['search_engine_score[1]']],
                        "evalue": _get_evalue_(psm, column_id)
                    },
                    "proteins":
                    proteins,
                    "num_missed_cleavages":
                    "-1"
                })
                if 'uri' in column_id.keys():
                    temp_hash['search_hit']['uri'] = psm[column_id['uri']]
            psm_hash.append(temp_hash)
    return psm_hash


#
# retrieve hit rank for a psm
#
def _get_hit_rank_(psm, column_id):
    '''
    :param psm: psm
    :param column_id:column number of rank
    :return: rank
    '''
    if 'rank' in column_id.keys():
        return psm[column_id['rank']]
    else:
        return "*"


#
# retrieve evalue from psm
#
def _get_evalue_(psm, column_id):
    '''
    :param psm: psm
    :param column_id: column number for evalue
    :return: evalue
    '''
    if 'fdr' in column_id.keys():
        return psm[column_id['fdr']]
    else:
        return '-1'


#
# calculates massdiff
#
def _calc_massdiff_(exp_mass, calc_mass):
    '''
    :param exp_mass: experimental mass
    :param calc_mass: calculated mass
    :return: mass difference
    '''
    if exp_mass == "" and calc_mass == "":
        return "*"
    else:
        return str(float(exp_mass) - float(calc_mass)) + ';' + str(
            exp_mass) + ';' + str(calc_mass)


#
# connect to unimod DB, parse and store in dict
#
def _unimod_parser_():
    '''
    :return: dictionairy with unimod IDs as keys and avg mass as value
    '''
    # create UNIMOD dictionairy (record_id->avg mass)
    unimod_dict = {}
    # connect UNIMOD database
    unimod_db = mass.Unimod(source='http://www.unimod.org/xml/unimod.xml')
    for mod in unimod_db.mods:
        unimod_dict[(mod['record_id'])] = mod['avge_mass']
    return unimod_dict


#
# Parse modifications (from UNIMOD/PSIMOD to neutral loss)
#
def _get_modifications_(mods):
    '''
    :param peptide: peptide sequence
    :param mods: peptide modification in unimod or psimod format
    :param unimod: unimod dictionairy
    :param psimod: psimod dictionairy
    :return: list of modification dictionairies for this peptide
    '''
    # seperate modification string into list of modifications
    mods = mods.replace(', ', ';')
    mod_list = mods.split(',')
    for i in range(0, len(mod_list)):
        mod_list[i].replace(';', ', ')

    # create variable to store results
    modification = []

    #iterate over all modifications
    for mod in mod_list:
        if mod == "0" or mod.upper() == 'NULL':
            break

        #partition unimod variable, separating position from unimod ID
        mod_partitions = mod.split('-', 1)

        #calculate mass
        if "UNIMOD" in mod_partitions[1].upper():
            #look up avg mass in unimod dict and store modification array
            modification.append({
                "position": mod_partitions[0],
                "mass": str(mod_partitions[1])
            })
        elif "MOD" in mod_partitions[1].upper():
            modification.append({
                "position": mod_partitions[0],
                "mass": str(mod_partitions[1])
            })

    return modification


#
# Parse modifications (from UNIMOD/PSIMOD to neutral loss)
#
def _get_modified_sequence_(peptide, mods, unimod, psimod):
    '''
    :param peptide: peptide sequence
    :param mods: peptide modification in unimod or psimod format
    :param unimod: unimod dictionairy
    :param psimod: psimod dictionairy
    :return: list of modification dictionairies for this peptide
    '''
    # seperate modification string into list of modifications
    mod_list = mods.split(',')
    mod_sequence = peptide

    # create variable to store results
    modification = []

    #iterate over all modifications
    for mod in mod_list:
        if mod == "0" or mod.upper() == 'NULL':
            break

        #partition unimod variable, separating position from unimod ID
        mod_partitions = mod.split("-")

        #calculate mass
        db_type = mod_partitions[1].split(':')[0]
        if db_type == "UNIMOD":
            #look up avg mass in unimod dict and store modification array
            mass = unimod[int(mod_partitions[1].split(':')[1])]
            modification.append([mod_partitions[0], mass])
        elif db_type == "MOD":
            mass = psimod[mod_partitions[1]]
            modification.append([mod_partitions[0], mass])
    modification = sorted(modification)
    shift = 0
    for mod in modification:
        pos = int(mod[0]) + shift + 1
        mod_bracketed = "[" + str(mod[1]) + "]"
        mod_sequence = mod_sequence[:pos] + mod_bracketed + mod_sequence[pos:]
        shift += len(mod_bracketed)
    return mod_sequence


#
# Parse modifications (from UNIMOD/PSIMOD to neutral loss)
#
def _get_modifications_neutral_(peptide, mods, unimod, psimod):
    '''
    :param peptide: peptide sequence
    :param mods: peptide modification in unimod or psimod format
    :param unimod: unimod dictionairy
    :param psimod: psimod dictionairy
    :return: list of modification dictionairies for this peptide
    '''
    # seperate modification string into list of modifications
    mod_list = mods.split(',')

    # create variable to store results
    modification = []

    #iterate over all modifications
    for mod in mod_list:
        if mod == "0" or mod.upper() == 'NULL':
            break

        #partition unimod variable, separating position from unimod ID
        mod_partitions = mod.split("-")

        #calculate mass
        db_type = mod_partitions[1].split(':')[0]
        if db_type == "UNIMOD":
            #look up avg mass in unimod dict and store modification array
            mass = unimod[int(mod_partitions[1].split(':')[1])]
            modification.append({"position": mod_partitions[0], "mass": mass})
        elif db_type == "MOD":
            mass = psimod[mod_partitions[1]]
            modification.append({"position": mod_partitions[0], "mass": mass})

    return modification


#
# extract comments from mzTAB
#
def extract_comments_from_mztab(psm_file):
    f = open(psm_file, 'r')
    comments = []
    for line in f:
        if line[0:3] == "COM" or line[0:3] == "MTD":
            comments.append(line)
        else:
            return comments
            break
    return comments


#
# Extract enzyme specificity from mzTab
#
def get_enzyme_specificity_mztab(psm_file):
    return 3


#
# Extract enzyme specificity from mzTab
#
def get_enzyme_mztab(psm_file):
    return "-1"


#
# Connect to psimod, parse IDs and store in dict
#
def _psimod_xml_parser_():
    '''
    :return: dictionairy with psimod IDs as key and avg mass as value
    '''
    # create dictionairy to store id/avg difference mass
    psimod = {}
    # parse and store psimod data
    import urllib2
    from lxml import etree
    data = urllib2.urlopen(
        "https://raw.githubusercontent.com/HUPO-PSI/psi-mod-CV/master/PSI-MOD.obo.xml"
    )
    root = etree.parse(data)
    root2 = etree.tostring(root)
    for term in root.xpath('//term'):
        for id in term.xpath('./id'):
            id = id.text
        for xref in term.xpath('./xref_analog'):
            for dbname in xref.xpath('./dbname'):
                if dbname.text == "DiffAvg":
                    avg_diff_mass = xref.xpath('./name/text()')
                    psimod[id] = avg_diff_mass[0]
    return psimod


#get_PSM_mztab("/home/vladie/Desktop/proBAMconvert/HCT116_NtermCofr.mztab")
