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

import proBAM_mzid
import proBAM_mzTab
import proBAM_pepxml
import re


#
# Import PSM file and parse into dictionairy
#
def get_PSM_hash(psm_file,decoy_annotation,validated_only):
    '''
    :param psm_file: psm file (pepxml,mzid or mztab)
    :raise: IO error: unrecognized file format
    :return: dictionairy of parsed psm file,
    '''
    print "Reading PSM-file"
    try:
        # catch mzid file format and parse
        if re.match('^.*\.(mzid)$',psm_file.lower())!=None:
            PSM=proBAM_mzid.get_PSM_mzid(psm_file,validated_only)

        # catch pepxml file format and parse
        elif re.match('^.*\.(pepxml|pep.xml|xml)$',psm_file.lower())!=None:
            PSM=proBAM_pepxml.get_PSM_pepxml(psm_file)


        # catch mztab file format and parse
        elif re.match('^.*\.(mztab)$',psm_file.lower())!=None:
            PSM=proBAM_mzTab.get_PSM_mztab(psm_file)

        else:
            raise IOError('Unrecognized file extension, \n ' \
                  'Accepted file extensions: .mzid/.pepxml/.pep.xml/.xml')

    except Exception as e:
        print "ERROR: Unable to  the PSM file : \n"
        print e.__doc__
        print e.message
        print "\nPlease confirm that the file is conform with the document specification." \
              "If the error keeps occuring contact the developers at https://github.com/Biobix/proBAMconvert/issues " \
              "and supply this error message along with the file and used settings"
        raise IOError()



    # Debug section
    '''
    count=0
    PSM=[]
    for row in PEP:
        print row
        if count==10000:
            break
        count+=1
        PSM.append(row)
    del PEP
    '''
    return PSM