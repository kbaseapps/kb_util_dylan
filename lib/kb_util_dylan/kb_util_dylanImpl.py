# -*- coding: utf-8 -*-
#BEGIN_HEADER
import os
import sys
import shutil
import hashlib
import subprocess
import requests
import re
import traceback
import uuid
from datetime import datetime
from pprint import pprint, pformat
import numpy as np
import gzip

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import generic_protein
from biokbase.workspace.client import Workspace as workspaceService
from requests_toolbelt import MultipartEncoder  # added
from biokbase.AbstractHandle.Client import AbstractHandle as HandleService  # added

# SDK Utils
from ReadsUtils.ReadsUtilsClient import ReadsUtils
from SetAPI.SetAPIClient import SetAPI
from KBaseReport.KBaseReportClient import KBaseReport

# silence whining
import requests
requests.packages.urllib3.disable_warnings()

#END_HEADER


class kb_util_dylan:
    '''
    Module Name:
    kb_util_dylan

    Module Description:
    ** A KBase module: kb_util_dylan
**
** This module contains basic utilities
    '''

    ######## WARNING FOR GEVENT USERS #######
    # Since asynchronous IO can lead to methods - even the same method -
    # interrupting each other, you must be *very* careful when using global
    # state. A method could easily clobber the state set by another while
    # the latter method is running.
    #########################################
    VERSION = "0.0.1"
    GIT_URL = "https://github.com/kbaseapps/kb_util_dylan.git"
    GIT_COMMIT_HASH = "268802f9ed7060a0136335df682a636c7b351b80"
    
    #BEGIN_CLASS_HEADER
    workspaceURL = None
    shockURL = None
    handleURL = None

    # target is a list for collecting log messages
    def log(self, target, message):
        # we should do something better here...
        if target is not None:
            target.append(message)
        print(message)
        sys.stdout.flush()

    def get_single_end_read_library(self, ws_data, ws_info, forward):
        pass

    def get_feature_set_seqs(self, ws_data, ws_info):
        pass

    def get_genome_feature_seqs(self, ws_data, ws_info):
        pass

    def get_genome_set_feature_seqs(self, ws_data, ws_info):
        pass

    # config contains contents of config file in a hash or None if it couldn't
    # be found
    def __init__(self, config):
        #BEGIN_CONSTRUCTOR
        self.workspaceURL = config['workspace-url']
        self.shockURL = config['shock-url']
        self.handleURL = config['handle-service-url']
        self.scratch = os.path.abspath(config['scratch'])
        # HACK!! temporary hack for issue where megahit fails on mac because of silent named pipe error
        #self.host_scratch = self.scratch
        self.scratch = os.path.join('/kb','module','local_scratch')
        # end hack
        if not os.path.exists(self.scratch):
            os.makedirs(self.scratch)

        #END_CONSTRUCTOR
        pass


    # Helper script borrowed from the transform service, logger removed
    #
    def upload_file_to_shock(self,
                             console,  # DEBUG
                             shock_service_url = None,
                             filePath = None,
                             ssl_verify = True,
                             token = None):
        """
        Use HTTP multi-part POST to save a file to a SHOCK instance.
        """
        self.log(console,"UPLOADING FILE "+filePath+" TO SHOCK")

        if token is None:
            raise Exception("Authentication token required!")

        #build the header
        header = dict()
        header["Authorization"] = "Oauth {0}".format(token)
        if filePath is None:
            raise Exception("No file given for upload to SHOCK!")

        dataFile = open(os.path.abspath(filePath), 'rb')
        m = MultipartEncoder(fields={'upload': (os.path.split(filePath)[-1], dataFile)})
        header['Content-Type'] = m.content_type

        #logger.info("Sending {0} to {1}".format(filePath,shock_service_url))
        try:
            response = requests.post(shock_service_url + "/node", headers=header, data=m, allow_redirects=True, verify=ssl_verify)
            dataFile.close()
        except:
            dataFile.close()
            raise
        if not response.ok:
            response.raise_for_status()
        result = response.json()
        if result['error']:
            raise Exception(result['error'][0])
        else:
            return result["data"]


    def upload_SingleEndLibrary_to_shock_and_ws (self,
                                                 ctx,
                                                 console,  # DEBUG
                                                 workspace_name,
                                                 obj_name,
                                                 file_path,
                                                 provenance,
                                                 sequencing_tech):

        self.log(console,'UPLOADING FILE '+file_path+' TO '+workspace_name+'/'+obj_name)

        # 1) upload files to shock
        token = ctx['token']
        forward_shock_file = self.upload_file_to_shock(
            console,  # DEBUG
            shock_service_url = self.shockURL,
            filePath = file_path,
            token = token
            )
        #pprint(forward_shock_file)
        self.log(console,'SHOCK UPLOAD DONE')

        # 2) create handle
        self.log(console,'GETTING HANDLE')
        hs = HandleService(url=self.handleURL, token=token)
        forward_handle = hs.persist_handle({
                                        'id' : forward_shock_file['id'], 
                                        'type' : 'shock',
                                        'url' : self.shockURL,
                                        'file_name': forward_shock_file['file']['name'],
                                        'remote_md5': forward_shock_file['file']['checksum']['md5']})

        
        # 3) save to WS
        self.log(console,'SAVING TO WORKSPACE')
        single_end_library = {
            'lib': {
                'file': {
                    'hid':forward_handle,
                    'file_name': forward_shock_file['file']['name'],
                    'id': forward_shock_file['id'],
                    'url': self.shockURL,
                    'type':'shock',
                    'remote_md5':forward_shock_file['file']['checksum']['md5']
                },
                'encoding':'UTF8',
                'type':'fasta',
                'size':forward_shock_file['file']['size']
            },
            'sequencing_tech':sequencing_tech
        }
        self.log(console,'GETTING WORKSPACE SERVICE OBJECT')
        ws = workspaceService(self.workspaceURL, token=ctx['token'])
        self.log(console,'SAVE OPERATION...')
        new_obj_info = ws.save_objects({
                        'workspace':workspace_name,
                        'objects':[
                            {
                                'type':'KBaseFile.SingleEndLibrary',
                                'data':single_end_library,
                                'name':obj_name,
                                'meta':{},
                                'provenance':provenance
                            }]
                        })
        self.log(console,'SAVED TO WORKSPACE')

        return new_obj_info[0]

    #END_CLASS_HEADER

    # config contains contents of config file in a hash or None if it couldn't
    # be found
    def __init__(self, config):
        #BEGIN_CONSTRUCTOR
        self.workspaceURL = config['workspace-url']
        self.shockURL = config['shock-url']
        self.handleURL = config['handle-service-url']
        self.scratch = os.path.abspath(config['scratch'])
        # HACK!! temporary hack for issue where megahit fails on mac because of silent named pipe error
        #self.host_scratch = self.scratch
        self.scratch = os.path.join('/kb','module','local_scratch')
        # end hack
        if not os.path.exists(self.scratch):
            os.makedirs(self.scratch)

        #END_CONSTRUCTOR
        pass
    

    def KButil_FASTQ_to_FASTA(self, ctx, params):
        """
        :param params: instance of type "KButil_FASTQ_to_FASTA_Params"
           (KButil_FASTQ_to_FASTA() ** ** Method for Converting a FASTQ
           SingleEndLibrary to a FASTA SingleEndLibrary) -> structure:
           parameter "workspace_name" of type "workspace_name" (** The
           workspace object refs are of form: ** **    objects =
           ws.get_objects([{'ref':
           params['workspace_id']+'/'+params['obj_name']}]) ** ** "ref" means
           the entire name combining the workspace id and the object name **
           "id" is a numerical identifier of the workspace or object, and
           should just be used for workspace ** "name" is a string identifier
           of a workspace or object.  This is received from Narrative.),
           parameter "input_name" of type "data_obj_name", parameter
           "output_name" of type "data_obj_name"
        :returns: instance of type "KButil_FASTQ_to_FASTA_Output" ->
           structure: parameter "report_name" of type "data_obj_name",
           parameter "report_ref" of type "data_obj_ref"
        """
        # ctx is the context object
        # return variables are: returnVal
        #BEGIN KButil_FASTQ_to_FASTA
        console = []
        invalid_msgs = []
        self.log(console,'Running KButil_FASTQ_to_FASTA with params=')
        self.log(console, "\n"+pformat(params))
        report = ''
#        report = 'Running KButil_FASTQ_to_FASTA with params='
#        report += "\n"+pformat(params)


        #### do some basic checks
        #
        if 'workspace_name' not in params:
            raise ValueError('workspace_name parameter is required')
        if 'input_name' not in params:
            raise ValueError('input_name parameter is required')
        if 'output_name' not in params:
            raise ValueError('output_name parameter is required')


        # Obtain the input object
        #
        forward_reads_file_compression = None
        sequencing_tech = 'N/A'
        try:
            ws = workspaceService(self.workspaceURL, token=ctx['token'])
            objects = ws.get_objects([{'ref': params['workspace_name']+'/'+params['input_name']}])
            data = objects[0]['data']
            info = objects[0]['info']
            # Object Info Contents
            # absolute ref = info[6] + '/' + info[0] + '/' + info[4]
            # 0 - obj_id objid
            # 1 - obj_name name
            # 2 - type_string type
            # 3 - timestamp save_date
            # 4 - int version
            # 5 - username saved_by
            # 6 - ws_id wsid
            # 7 - ws_name workspace
            # 8 - string chsum
            # 9 - int size 
            # 10 - usermeta meta
            type_name = info[2].split('.')[1].split('-')[0]

            if type_name == 'SingleEndLibrary':
                type_namespace = info[2].split('.')[0]
                if type_namespace == 'KBaseAssembly':
                    file_name = data['handle']['file_name']
                elif type_namespace == 'KBaseFile':
                    file_name = data['lib']['file']['file_name']
                else:
                    raise ValueError('bad data type namespace: '+type_namespace)
                #self.log(console, 'INPUT_FILENAME: '+file_name)  # DEBUG
                if file_name[-3:] == ".gz":
                    forward_reads_file_compression = 'gz'
                if 'sequencing_tech' in data:
                    sequencing_tech = data['sequencing_tech']

        except Exception as e:
            raise ValueError('Unable to fetch input_name object from workspace: ' + str(e))
            #to get the full stack trace: traceback.format_exc()
        
        # pull data from SHOCK
        #
        try:
            if 'lib' in data:
                forward_reads = data['lib']['file']
            elif 'handle' in data:
                forward_reads = data['handle']
            else:
                self.log(console,"bad structure for 'forward_reads'")
                raise ValueError("bad structure for 'forward_reads'")

            ### NOTE: this section is what could be replaced by the transform services
            forward_reads_file_path = os.path.join(self.scratch,forward_reads['file_name'])
            forward_reads_file_handle = open(forward_reads_file_path, 'w', 0)
            self.log(console, 'downloading reads file: '+str(forward_reads_file_path))
            headers = {'Authorization': 'OAuth '+ctx['token']}
            r = requests.get(forward_reads['url']+'/node/'+forward_reads['id']+'?download', stream=True, headers=headers)
            for chunk in r.iter_content(1024):
                forward_reads_file_handle.write(chunk)
            forward_reads_file_handle.close();
            self.log(console, 'done')
            ### END NOTE
        except Exception as e:
            print(traceback.format_exc())
            raise ValueError('Unable to download single-end read library files: ' + str(e))


        #### Create the file to upload
        ##
        output_file_name   = params['output_name']+'.fna'
        output_file_path  = os.path.join(self.scratch,output_file_name)
        input_file_handle  = open(forward_reads_file_path, 'r', -1)
        output_file_handle = open(output_file_path, 'w', -1)
        self.log(console, 'PROCESSING reads file: '+str(forward_reads_file_path))

        seq_cnt = 0
        header = None
        last_header = None
        last_seq_buf = None
        last_line_was_header = False
        for line in input_file_handle:
            if line.startswith('@'):
                seq_cnt += 1
                header = line[1:]
                if last_header != None:
                    output_file_handle.write('>'+last_header)
                    output_file_handle.write(last_seq_buf)
                last_seq_buf = None
                last_header = header
                last_line_was_header = True
            elif last_line_was_header:
                last_seq_buf = line
                last_line_was_header = False
            else:
                continue
        if last_header != None:
            output_file_handle.write('>'+last_header)
            output_file_handle.write(last_seq_buf)

        input_file_handle.close()
        output_file_handle.close()
        

        # load the method provenance from the context object
        #
        self.log(console,"SETTING PROVENANCE")  # DEBUG
        provenance = [{}]
        if 'provenance' in ctx:
            provenance = ctx['provenance']
        # add additional info to provenance here, in this case the input data object reference
        provenance[0]['input_ws_objects'] = []
        provenance[0]['input_ws_objects'].append(params['workspace_name']+'/'+params['input_name'])
        provenance[0]['service'] = 'kb_util_dylan'
        provenance[0]['method'] = 'KButil_FASTQ_to_FASTA'


        # Upload results
        #
        if len(invalid_msgs) == 0:
            self.log(console,"UPLOADING RESULTS")  # DEBUG
            self.upload_SingleEndLibrary_to_shock_and_ws (ctx,
                                                      console,  # DEBUG
                                                      params['workspace_name'],
                                                      params['output_name'],
                                                      output_file_path,
                                                      provenance,
                                                      sequencing_tech
                                                      )

        # build output report object
        #
        self.log(console,"BUILDING REPORT")  # DEBUG
        if len(invalid_msgs) == 0:
            report += 'sequences in library:  '+str(seq_cnt)+"\n"
            reportObj = {
                'objects_created':[{'ref':params['workspace_name']+'/'+params['output_name'], 'description':'KButil_FASTQ_to_FASTA'}],
                'text_message':report
                }
        else:
            report += "FAILURE:\n\n"+"\n".join(invalid_msgs)+"\n"
            reportObj = {
                'objects_created':[],
                'text_message':report
                }

        reportName = 'kbutil_fastq_to_fasta_report_'+str(uuid.uuid4())
        ws = workspaceService(self.workspaceURL, token=ctx['token'])
        report_obj_info = ws.save_objects({
#                'id':info[6],
                'workspace':params['workspace_name'],
                'objects':[
                    {
                        'type':'KBaseReport.Report',
                        'data':reportObj,
                        'name':reportName,
                        'meta':{},
                        'hidden':1,
                        'provenance':provenance
                    }
                ]
            })[0]

        self.log(console,"BUILDING RETURN OBJECT")
#        returnVal = { 'output_report_name': reportName,
#                      'output_report_ref': str(report_obj_info[6]) + '/' + str(report_obj_info[0]) + '/' + str(report_obj_info[4]),
#                      'output_filtered_ref': params['workspace_name']+'/'+params['output_filtered_name']
#                      }
        returnVal = { 'report_name': reportName,
                      'report_ref': str(report_obj_info[6]) + '/' + str(report_obj_info[0]) + '/' + str(report_obj_info[4]),
                      }
        self.log(console,"KButil_FASTQ_to_FASTA DONE")

        #END KButil_FASTQ_to_FASTA

        # At some point might do deeper type checking...
        if not isinstance(returnVal, dict):
            raise ValueError('Method KButil_FASTQ_to_FASTA return value ' +
                             'returnVal is not type dict as required.')
        # return the results
        return [returnVal]

    def KButil_Merge_FeatureSet_Collection(self, ctx, params):
        """
        :param params: instance of type
           "KButil_Merge_FeatureSet_Collection_Params"
           (KButil_Merge_FeatureSet_Collection() ** **  Method for merging
           FeatureSets) -> structure: parameter "workspace_name" of type
           "workspace_name" (** The workspace object refs are of form: ** ** 
           objects = ws.get_objects([{'ref':
           params['workspace_id']+'/'+params['obj_name']}]) ** ** "ref" means
           the entire name combining the workspace id and the object name **
           "id" is a numerical identifier of the workspace or object, and
           should just be used for workspace ** "name" is a string identifier
           of a workspace or object.  This is received from Narrative.),
           parameter "input_names" of type "data_obj_name", parameter
           "output_name" of type "data_obj_name", parameter "desc" of String
        :returns: instance of type
           "KButil_Merge_FeatureSet_Collection_Output" -> structure:
           parameter "report_name" of type "data_obj_name", parameter
           "report_ref" of type "data_obj_ref"
        """
        # ctx is the context object
        # return variables are: returnVal
        #BEGIN KButil_Merge_FeatureSet_Collection
        console = []
        invalid_msgs = []
        self.log(console,'Running KButil_Merge_FeatureSet_Collection with params=')
        self.log(console, "\n"+pformat(params))
        report = ''
#        report = 'Running KButil_FASTQ_to_FASTA with params='
#        report += "\n"+pformat(params)


        #### do some basic checks
        #
        if 'workspace_name' not in params:
            raise ValueError('workspace_name parameter is required')
        if 'desc' not in params:
            raise ValueError('desc parameter is required')
        if 'input_names' not in params:
            raise ValueError('input_names parameter is required')
        if 'output_name' not in params:
            raise ValueError('output_name parameter is required')


        # Build FeatureSet
        #
        element_ordering = []
        elements = {}
        featureSet_seen = dict()
        for featureSet_name in params['input_names']:
            if not featureSet_name in featureSet_seen.keys():
                featureSet_seen[featureSet_name] = 1
            else:
                self.log("repeat featureSet_name: '"+featureSet_name+"'")
                self.log(invalid_msgs,"repeat featureSet_name: '"+featureSet_name+"'")
                continue

            try:
                ws = workspaceService(self.workspaceURL, token=ctx['token'])
                objects = ws.get_objects([{'ref': params['workspace_name']+'/'+featureSet_name}])
                data = objects[0]['data']
                info = objects[0]['info']
                # Object Info Contents
                # absolute ref = info[6] + '/' + info[0] + '/' + info[4]
                # 0 - obj_id objid
                # 1 - obj_name name
                # 2 - type_string type
                # 3 - timestamp save_date
                # 4 - int version
                # 5 - username saved_by
                # 6 - ws_id wsid
                # 7 - ws_name workspace
                # 8 - string chsum
                # 9 - int size 
                # 10 - usermeta meta
                type_name = info[2].split('.')[1].split('-')[0]

            except Exception as e:
                raise ValueError('Unable to fetch input_name object from workspace: ' + str(e))
                #to get the full stack trace: traceback.format_exc()

            if type_name != 'FeatureSet':
                raise ValueError("Bad Type:  Should be FeatureSet instead of '"+type_name+"'")

            this_featureSet = data
            this_element_ordering = []
            if 'element_ordering' in this_featureSet.keys():
                this_element_ordering = this_featureSet['element_ordering']
            else:
                this_element_ordering = sorted(this_featureSet['elements'].keys())
            element_ordering.extend(this_element_ordering)
            self.log(console,'features in input set '+featureSet_name+': '+str(len(this_element_ordering)))
            report += 'features in input set '+featureSet_name+': '+str(len(this_element_ordering))+"\n"

            for fId in this_featureSet['elements'].keys():
                elements[fId] = this_featureSet['elements'][fId]
            

        # load the method provenance from the context object
        #
        self.log(console,"SETTING PROVENANCE")  # DEBUG
        provenance = [{}]
        if 'provenance' in ctx:
            provenance = ctx['provenance']
        # add additional info to provenance here, in this case the input data object reference
        provenance[0]['input_ws_objects'] = []
        for featureSet_name in params['input_names']:
            provenance[0]['input_ws_objects'].append(params['workspace_name']+'/'+featureSet_name)
        provenance[0]['service'] = 'kb_util_dylan'
        provenance[0]['method'] = 'KButil_Merge_FeatureSet_Collection'


        # Store output object
        #
        if len(invalid_msgs) == 0:
            self.log(console,"SAVING FEATURESET")  # DEBUG
            output_FeatureSet = {
                              'description': params['desc'],
                              'element_ordering': element_ordering,
                              'elements': elements
                            }

            new_obj_info = ws.save_objects({
                            'workspace': params['workspace_name'],
                            'objects':[{
                                    'type': 'KBaseCollections.FeatureSet',
                                    'data': output_FeatureSet,
                                    'name': params['output_name'],
                                    'meta': {},
                                    'provenance': provenance
                                }]
                        })


        # build output report object
        #
        self.log(console,"BUILDING REPORT")  # DEBUG
        if len(invalid_msgs) == 0:
            self.log(console,"features in output set "+params['output_name']+": "+str(len(element_ordering)))
            report += 'features in output set '+params['output_name']+': '+str(len(element_ordering))+"\n"
            reportObj = {
                'objects_created':[{'ref':params['workspace_name']+'/'+params['output_name'], 'description':'KButil_Merge_FeatureSet_Collection'}],
                'text_message':report
                }
        else:
            report += "FAILURE:\n\n"+"\n".join(invalid_msgs)+"\n"
            reportObj = {
                'objects_created':[],
                'text_message':report
                }

        reportName = 'kb_util_dylan_merge_featureset_report_'+str(uuid.uuid4())
        ws = workspaceService(self.workspaceURL, token=ctx['token'])
        report_obj_info = ws.save_objects({
#                'id':info[6],
                'workspace':params['workspace_name'],
                'objects':[
                    {
                        'type':'KBaseReport.Report',
                        'data':reportObj,
                        'name':reportName,
                        'meta':{},
                        'hidden':1,
                        'provenance':provenance
                    }
                ]
            })[0]


        # Build report and return
        #
        self.log(console,"BUILDING RETURN OBJECT")
        returnVal = { 'report_name': reportName,
                      'report_ref': str(report_obj_info[6]) + '/' + str(report_obj_info[0]) + '/' + str(report_obj_info[4]),
                      }
        self.log(console,"KButil_Merge_FeatureSet_Collection DONE")
        #END KButil_Merge_FeatureSet_Collection

        # At some point might do deeper type checking...
        if not isinstance(returnVal, dict):
            raise ValueError('Method KButil_Merge_FeatureSet_Collection return value ' +
                             'returnVal is not type dict as required.')
        # return the results
        return [returnVal]

    def KButil_Merge_GenomeSets(self, ctx, params):
        """
        :param params: instance of type "KButil_Merge_GenomeSets_Params"
           (KButil_Merge_GenomeSets() ** **  Method for merging GenomeSets)
           -> structure: parameter "workspace_name" of type "workspace_name"
           (** The workspace object refs are of form: ** **    objects =
           ws.get_objects([{'ref':
           params['workspace_id']+'/'+params['obj_name']}]) ** ** "ref" means
           the entire name combining the workspace id and the object name **
           "id" is a numerical identifier of the workspace or object, and
           should just be used for workspace ** "name" is a string identifier
           of a workspace or object.  This is received from Narrative.),
           parameter "input_names" of type "data_obj_name", parameter
           "output_name" of type "data_obj_name", parameter "desc" of String
        :returns: instance of type "KButil_Merge_GenomeSets_Output" ->
           structure: parameter "report_name" of type "data_obj_name",
           parameter "report_ref" of type "data_obj_ref"
        """
        # ctx is the context object
        # return variables are: returnVal
        #BEGIN KButil_Merge_GenomeSets
        console = []
        invalid_msgs = []
        self.log(console,'Running KButil_Merge_GenomeSets with params=')
        self.log(console, "\n"+pformat(params))
        report = ''
#        report = 'Running KButil_Merge_GenomeSets with params='
#        report += "\n"+pformat(params)


        #### do some basic checks
        #
        if 'workspace_name' not in params:
            raise ValueError('workspace_name parameter is required')
        if 'desc' not in params:
            raise ValueError('desc parameter is required')
        if 'input_names' not in params:
            raise ValueError('input_names parameter is required')
        if 'output_name' not in params:
            raise ValueError('output_name parameter is required')


        # Build GenomeSet
        #
        elements = dict()


        # Add Genomes from GenomeSets
        #
        for input_genomeset_name in params['input_names']:

            try:
                ws = workspaceService(self.workspaceURL, token=ctx['token'])
                objects = ws.get_objects([{'ref': params['workspace_name']+'/'+input_genomeset_name}])
                genomeSet = objects[0]['data']
                info = objects[0]['info']
                
                type_name = info[2].split('.')[1].split('-')[0]
                if type_name != 'GenomeSet':
                    raise ValueError("Bad Type:  Should be GenomeSet instead of '"+type_name+"'")
            except Exception as e:
                raise ValueError('Unable to fetch input_name object from workspace: ' + str(e))
                #to get the full stack trace: traceback.format_exc()

            for gId in genomeSet['elements'].keys():
                genomeRef = genomeSet['elements'][gId]['ref']
                try:
                    already_included = elements[gId]
                except:
                    elements[gId] = dict()
                    elements[gId]['ref'] = genomeRef  # the key line
                    self.log(console,"adding element "+gId+" : "+genomeRef)  # DEBUG
            

        # load the method provenance from the context object
        #
        self.log(console,"SETTING PROVENANCE")  # DEBUG
        provenance = [{}]
        if 'provenance' in ctx:
            provenance = ctx['provenance']
        # add additional info to provenance here, in this case the input data object reference
        try:
            prov_defined = provenance[0]['input_ws_objects']
        except:
            provenance[0]['input_ws_objects'] = []
        provenance[0]['input_ws_objects'].append(params['workspace_name']+'/'+params['input_genome_names'])
        provenance[0]['input_ws_objects'].append(params['workspace_name']+'/'+params['input_genomeset_name'])
        provenance[0]['service'] = 'kb_util_dylan'
        provenance[0]['method'] = 'KButil_Merge_GenomeSets'


        # Store output object
        #
        if len(invalid_msgs) == 0:
            self.log(console,"SAVING GENOMESET")  # DEBUG
            output_GenomeSet = {
                              'description': params['desc'],
                              'elements': elements
                            }

            new_obj_info = ws.save_objects({
                            'workspace': params['workspace_name'],
                            'objects':[{
                                    'type': 'KBaseSearch.GenomeSet',
                                    'data': output_GenomeSet,
                                    'name': params['output_name'],
                                    'meta': {},
                                    'provenance': provenance
                                }]
                        })


        # build output report object
        #
        self.log(console,"BUILDING REPORT")  # DEBUG
        if len(invalid_msgs) == 0:
            self.log(console,"genomes in output set "+params['output_name']+": "+str(len(elements.keys())))
            report += 'genomes in output set '+params['output_name']+': '+str(len(elements.keys()))+"\n"
            reportObj = {
                'objects_created':[{'ref':params['workspace_name']+'/'+params['output_name'], 'description':'KButil_Merge_GenomeSets'}],
                'text_message':report
                }
        else:
            report += "FAILURE:\n\n"+"\n".join(invalid_msgs)+"\n"
            reportObj = {
                'objects_created':[],
                'text_message':report
                }

        reportName = 'kb_util_dylan_merge_genomesets_report_'+str(uuid.uuid4())
        ws = workspaceService(self.workspaceURL, token=ctx['token'])
        report_obj_info = ws.save_objects({
#                'id':info[6],
                'workspace':params['workspace_name'],
                'objects':[
                    {
                        'type':'KBaseReport.Report',
                        'data':reportObj,
                        'name':reportName,
                        'meta':{},
                        'hidden':1,
                        'provenance':provenance
                    }
                ]
            })[0]


        # Build report and return
        #
        self.log(console,"BUILDING RETURN OBJECT")
        returnVal = { 'report_name': reportName,
                      'report_ref': str(report_obj_info[6]) + '/' + str(report_obj_info[0]) + '/' + str(report_obj_info[4]),
                      }
        self.log(console,"KButil_Merge_GenomeSets DONE")
        #END KButil_Merge_GenomeSets

        # At some point might do deeper type checking...
        if not isinstance(returnVal, dict):
            raise ValueError('Method KButil_Merge_GenomeSets return value ' +
                             'returnVal is not type dict as required.')
        # return the results
        return [returnVal]

    def KButil_Build_GenomeSet(self, ctx, params):
        """
        :param params: instance of type "KButil_Build_GenomeSet_Params"
           (KButil_Build_GenomeSet() ** **  Method for creating a GenomeSet)
           -> structure: parameter "workspace_name" of type "workspace_name"
           (** The workspace object refs are of form: ** **    objects =
           ws.get_objects([{'ref':
           params['workspace_id']+'/'+params['obj_name']}]) ** ** "ref" means
           the entire name combining the workspace id and the object name **
           "id" is a numerical identifier of the workspace or object, and
           should just be used for workspace ** "name" is a string identifier
           of a workspace or object.  This is received from Narrative.),
           parameter "input_names" of type "data_obj_name", parameter
           "output_name" of type "data_obj_name", parameter "desc" of String
        :returns: instance of type "KButil_Build_GenomeSet_Output" ->
           structure: parameter "report_name" of type "data_obj_name",
           parameter "report_ref" of type "data_obj_ref"
        """
        # ctx is the context object
        # return variables are: returnVal
        #BEGIN KButil_Build_GenomeSet
        console = []
        invalid_msgs = []
        self.log(console,'Running KButil_Build_GenomeSet with params=')
        self.log(console, "\n"+pformat(params))
        report = ''
#        report = 'Running KButil_Build_GenomeSet with params='
#        report += "\n"+pformat(params)


        #### do some basic checks
        #
        if 'workspace_name' not in params:
            raise ValueError('workspace_name parameter is required')
        if 'desc' not in params:
            raise ValueError('desc parameter is required')
        if 'input_names' not in params:
            raise ValueError('input_names parameter is required')
        if 'output_name' not in params:
            raise ValueError('output_name parameter is required')


        # Build GenomeSet
        #
        elements = {}
        genome_seen = dict()
        
        for genome_name in params['input_names']:
            genomeRef = params['workspace_name'] + '/' + genome_name

            try:
                already_included = genome_seen[genomeRef]
            except:
                genome_seen[genomeRef] = True

                try:
                    ws = workspaceService(self.workspaceURL, token=ctx['token'])
                    objects = ws.get_objects([{'ref': genomeRef}])
                    data = objects[0]['data']
                    info = objects[0]['info']
                    genomeObj = data
                    type_name = info[2].split('.')[1].split('-')[0]
                except Exception as e:
                    raise ValueError('Unable to fetch input_name object from workspace: ' + str(e))
                if type_name != 'Genome' and type_name != 'GenomeAnnotation':
                    raise ValueError("Bad Type:  Should be Genome or GenomeAnnotation instead of '"+type_name+"' for ref: '"+genomeRef+"'")
                    
                genome_id = genomeObj['id'] if type_name == 'Genome' else genomeObj['genome_annotation_id']
                if not genome_id in elements.keys():
                    elements[genome_id] = dict()
                elements[genome_id]['ref'] = genomeRef  # the key line
                self.log(console,"adding element "+genome_id+" : "+genomeRef)  # DEBUG


        # load the method provenance from the context object
        #
        self.log(console,"SETTING PROVENANCE")  # DEBUG
        provenance = [{}]
        if 'provenance' in ctx:
            provenance = ctx['provenance']
        # add additional info to provenance here, in this case the input data object reference
        provenance[0]['input_ws_objects'] = []
        for genome_name in params['input_names']:
            provenance[0]['input_ws_objects'].append(params['workspace_name']+'/'+genome_name)
        provenance[0]['service'] = 'kb_util_dylan'
        provenance[0]['method'] = 'KButil_Build_GenomeSet'


        # Store output object
        #
        if len(invalid_msgs) == 0:
            self.log(console,"SAVING GENOMESET")  # DEBUG
            output_GenomeSet = {
                              'description': params['desc'],
                              'elements': elements
                            }

            new_obj_info = ws.save_objects({
                            'workspace': params['workspace_name'],
                            'objects':[{
                                    'type': 'KBaseSearch.GenomeSet',
                                    'data': output_GenomeSet,
                                    'name': params['output_name'],
                                    'meta': {},
                                    'provenance': provenance
                                }]
                        })


        # build output report object
        #
        self.log(console,"BUILDING REPORT")  # DEBUG
        if len(invalid_msgs) == 0:
            self.log(console,"genomes in output set "+params['output_name']+": "+str(len(elements.keys())))
            report += 'genomes in output set '+params['output_name']+': '+str(len(elements.keys()))+"\n"
            reportObj = {
                'objects_created':[{'ref':params['workspace_name']+'/'+params['output_name'], 'description':'KButil_Build_GenomeSet'}],
                'text_message':report
                }
        else:
            report += "FAILURE:\n\n"+"\n".join(invalid_msgs)+"\n"
            reportObj = {
                'objects_created':[],
                'text_message':report
                }

        reportName = 'kb_util_dylan_build_genomeset_report_'+str(uuid.uuid4())
        ws = workspaceService(self.workspaceURL, token=ctx['token'])
        report_obj_info = ws.save_objects({
#                'id':info[6],
                'workspace':params['workspace_name'],
                'objects':[
                    {
                        'type':'KBaseReport.Report',
                        'data':reportObj,
                        'name':reportName,
                        'meta':{},
                        'hidden':1,
                        'provenance':provenance
                    }
                ]
            })[0]


        # Build report and return
        #
        self.log(console,"BUILDING RETURN OBJECT")
        returnVal = { 'report_name': reportName,
                      'report_ref': str(report_obj_info[6]) + '/' + str(report_obj_info[0]) + '/' + str(report_obj_info[4]),
                      }
        self.log(console,"KButil_Build_GenomeSet DONE")
        #END KButil_Build_GenomeSet

        # At some point might do deeper type checking...
        if not isinstance(returnVal, dict):
            raise ValueError('Method KButil_Build_GenomeSet return value ' +
                             'returnVal is not type dict as required.')
        # return the results
        return [returnVal]

    def KButil_Build_GenomeSet_from_FeatureSet(self, ctx, params):
        """
        :param params: instance of type
           "KButil_Build_GenomeSet_from_FeatureSet_Params"
           (KButil_Build_GenomeSet_from_FeatureSet() ** **  Method for
           obtaining a GenomeSet from a FeatureSet) -> structure: parameter
           "workspace_name" of type "workspace_name" (** The workspace object
           refs are of form: ** **    objects = ws.get_objects([{'ref':
           params['workspace_id']+'/'+params['obj_name']}]) ** ** "ref" means
           the entire name combining the workspace id and the object name **
           "id" is a numerical identifier of the workspace or object, and
           should just be used for workspace ** "name" is a string identifier
           of a workspace or object.  This is received from Narrative.),
           parameter "input_name" of type "data_obj_name", parameter
           "output_name" of type "data_obj_name", parameter "desc" of String
        :returns: instance of type
           "KButil_Build_GenomeSet_from_FeatureSet_Output" -> structure:
           parameter "report_name" of type "data_obj_name", parameter
           "report_ref" of type "data_obj_ref"
        """
        # ctx is the context object
        # return variables are: returnVal
        #BEGIN KButil_Build_GenomeSet_from_FeatureSet
        console = []
        invalid_msgs = []
        self.log(console,'Running KButil_Build_GenomeSet_from_FeatureSet with params=')
        self.log(console, "\n"+pformat(params))
        report = ''
#        report = 'Running KButil_Build_GenomeSet_from_FeatureSet with params='
#        report += "\n"+pformat(params)


        #### do some basic checks
        #
        if 'workspace_name' not in params:
            raise ValueError('workspace_name parameter is required')
        if 'desc' not in params:
            raise ValueError('desc parameter is required')
        if 'input_name' not in params:
            raise ValueError('input_name parameter is required')
        if 'output_name' not in params:
            raise ValueError('output_name parameter is required')


        # Obtain FeatureSet
        #
        try:
            ws = workspaceService(self.workspaceURL, token=ctx['token'])
            objects = ws.get_objects([{'ref': params['workspace_name']+'/'+params['input_name']}])
            data = objects[0]['data']
            info = objects[0]['info']
            # Object Info Contents
            # absolute ref = info[6] + '/' + info[0] + '/' + info[4]
            # 0 - obj_id objid
            # 1 - obj_name name
            # 2 - type_string type
            # 3 - timestamp save_date
            # 4 - int version
            # 5 - username saved_by
            # 6 - ws_id wsid
            # 7 - ws_name workspace
            # 8 - string chsum
            # 9 - int size 
            # 10 - usermeta meta
            featureSet = data
            type_name = info[2].split('.')[1].split('-')[0]
        except Exception as e:
            raise ValueError('Unable to fetch input_name object from workspace: ' + str(e))
            #to get the full stack trace: traceback.format_exc()
        if type_name != 'FeatureSet':
            raise ValueError("Bad Type:  Should be FeatureSet instead of '"+type_name+"'")


        # Build GenomeSet
        #
        elements = {}
        genome_seen = dict()
        
        for fId in featureSet['elements'].keys():
            for genomeRef in featureSet['elements'][fId]:

                try:
                    already_included = genome_seen[genomeRef]
                except:
                    genome_seen[genomeRef] = True

                    try:
                        ws = workspaceService(self.workspaceURL, token=ctx['token'])
                        objects = ws.get_objects([{'ref': genomeRef}])
                        data = objects[0]['data']
                        info = objects[0]['info']
                        genomeObj = data
                        type_name = info[2].split('.')[1].split('-')[0]
                    except Exception as e:
                        raise ValueError('Unable to fetch input_name object from workspace: ' + str(e))
                    if type_name != 'Genome' and type_name != 'GenomeAnnotaton':
                        raise ValueError("Bad Type:  Should be Genome or GenomeAnnotation instead of '"+type_name+"' for ref: '"+genomeRef+"'")
                    
                    genome_id = genomeObj['id'] if type_name == 'Genome' else genomeObj['genome_annotation_id']
                    if not genome_id in elements.keys():
                        elements[genome_id] = dict()
                    elements[genome_id]['ref'] = genomeRef  # the key line
                    self.log(console,"adding element "+genome_id+" : "+genomeRef)  # DEBUG
            

        # load the method provenance from the context object
        #
        self.log(console,"SETTING PROVENANCE")  # DEBUG
        provenance = [{}]
        if 'provenance' in ctx:
            provenance = ctx['provenance']
        # add additional info to provenance here, in this case the input data object reference
        provenance[0]['input_ws_objects'] = []
        provenance[0]['input_ws_objects'].append(params['workspace_name']+'/'+params['input_name'])
        provenance[0]['service'] = 'kb_util_dylan'
        provenance[0]['method'] = 'KButil_Build_GenomeSet_from_FeatureSet'


        # Store output object
        #
        if len(invalid_msgs) == 0:
            self.log(console,"SAVING GENOMESET")  # DEBUG
            output_GenomeSet = {
                              'description': params['desc'],
                              'elements': elements
                            }

            new_obj_info = ws.save_objects({
                            'workspace': params['workspace_name'],
                            'objects':[{
                                    'type': 'KBaseSearch.GenomeSet',
                                    'data': output_GenomeSet,
                                    'name': params['output_name'],
                                    'meta': {},
                                    'provenance': provenance
                                }]
                        })


        # build output report object
        #
        self.log(console,"BUILDING REPORT")  # DEBUG
        if len(invalid_msgs) == 0:
            self.log(console,"genomes in output set "+params['output_name']+": "+str(len(elements.keys())))
            report += 'genomes in output set '+params['output_name']+': '+str(len(elements.keys()))+"\n"
            reportObj = {
                'objects_created':[{'ref':params['workspace_name']+'/'+params['output_name'], 'description':'KButil_Build_GenomeSet_from_FeatureSet'}],
                'text_message':report
                }
        else:
            report += "FAILURE:\n\n"+"\n".join(invalid_msgs)+"\n"
            reportObj = {
                'objects_created':[],
                'text_message':report
                }

        reportName = 'kb_util_dylan_build_genomeset_from_featureset_report_'+str(uuid.uuid4())
        ws = workspaceService(self.workspaceURL, token=ctx['token'])
        report_obj_info = ws.save_objects({
#                'id':info[6],
                'workspace':params['workspace_name'],
                'objects':[
                    {
                        'type':'KBaseReport.Report',
                        'data':reportObj,
                        'name':reportName,
                        'meta':{},
                        'hidden':1,
                        'provenance':provenance
                    }
                ]
            })[0]


        # Build report and return
        #
        self.log(console,"BUILDING RETURN OBJECT")
        returnVal = { 'report_name': reportName,
                      'report_ref': str(report_obj_info[6]) + '/' + str(report_obj_info[0]) + '/' + str(report_obj_info[4]),
                      }
        self.log(console,"KButil_Build_GenomeSet_from_FeatureSet DONE")
        #END KButil_Build_GenomeSet_from_FeatureSet

        # At some point might do deeper type checking...
        if not isinstance(returnVal, dict):
            raise ValueError('Method KButil_Build_GenomeSet_from_FeatureSet return value ' +
                             'returnVal is not type dict as required.')
        # return the results
        return [returnVal]

    def KButil_Add_Genomes_to_GenomeSet(self, ctx, params):
        """
        :param params: instance of type
           "KButil_Add_Genomes_to_GenomeSet_Params"
           (KButil_Add_Genomes_to_GenomeSet() ** **  Method for adding a
           Genome to a GenomeSet) -> structure: parameter "workspace_name" of
           type "workspace_name" (** The workspace object refs are of form:
           ** **    objects = ws.get_objects([{'ref':
           params['workspace_id']+'/'+params['obj_name']}]) ** ** "ref" means
           the entire name combining the workspace id and the object name **
           "id" is a numerical identifier of the workspace or object, and
           should just be used for workspace ** "name" is a string identifier
           of a workspace or object.  This is received from Narrative.),
           parameter "input_genome_names" of type "data_obj_name", parameter
           "input_genomeset_name" of type "data_obj_name", parameter
           "output_name" of type "data_obj_name", parameter "desc" of String
        :returns: instance of type "KButil_Add_Genomes_to_GenomeSet_Output"
           -> structure: parameter "report_name" of type "data_obj_name",
           parameter "report_ref" of type "data_obj_ref"
        """
        # ctx is the context object
        # return variables are: returnVal
        #BEGIN KButil_Add_Genomes_to_GenomeSet
        console = []
        invalid_msgs = []
        self.log(console,'Running KButil_Add_Genomes_to_GenomeSet with params=')
        self.log(console, "\n"+pformat(params))
        report = ''
#        report = 'Running KButil_Add_Genomes_to_GenomeSet with params='
#        report += "\n"+pformat(params)


        #### do some basic checks
        #
        if 'workspace_name' not in params:
            raise ValueError('workspace_name parameter is required')
        if 'desc' not in params:
            raise ValueError('desc parameter is required')
        if 'input_genome_names' not in params:
            raise ValueError('input_genome_names parameter is required')
        if 'input_genomeset_name' not in params:
            raise ValueError('input_genomeset_name parameter is required')
        if 'output_name' not in params:
            raise ValueError('output_name parameter is required')


        # Build GenomeSet
        #
        elements = dict()
        genome_seen = dict()

        # add new genome
        for genome_name in params['input_genome_names']:

            try:
                ws = workspaceService(self.workspaceURL, token=ctx['token'])
                objects = ws.get_objects([{'ref': params['workspace_name']+'/'+params['genome_name']}])
                genomeObj = objects[0]['data']
                info = objects[0]['info']

                genomeRef = str(info[6]) + '/' + str(info[0]) + '/' + str(info[4])
                type_name = info[2].split('.')[1].split('-')[0]
                if type_name != 'Genome' and type != 'GenomeAnnotation':
                    raise ValueError("Bad Type:  Should be Genome or GenomeAnnotation instead of '"+type_name+"'")

            except Exception as e:
                raise ValueError('Unable to fetch input_name object from workspace: ' + str(e))
                #to get the full stack trace: traceback.format_exc()
            
            gId = genomeObj['id'] if type_name == 'Genome' else genomeObj['genome_annotation_id']
            try:
                already_included = elements[gId]
            except:
                elements[gId] = dict()
                elements[gId]['ref'] = genomeRef  # the key line
                self.log(console,"adding new element "+gId+" : "+genomeRef)  # DEBUG


        # add rest of old GenomeSet
        #
        if 'input_genomeset_name' in params and params['input_genomeset_name'] != None:
            try:
                ws = workspaceService(self.workspaceURL, token=ctx['token'])
                objects = ws.get_objects([{'ref': params['workspace_name']+'/'+params['input_genomeset_name']}])
                genomeSet = objects[0]['data']
                info = objects[0]['info']
                
                type_name = info[2].split('.')[1].split('-')[0]
                if type_name != 'GenomeSet':
                    raise ValueError("Bad Type:  Should be GenomeSet instead of '"+type_name+"'")
            except Exception as e:
                raise ValueError('Unable to fetch input_name object from workspace: ' + str(e))
                #to get the full stack trace: traceback.format_exc()

            for gId in genomeSet['elements'].keys():
                genomeRef = genomeSet['elements'][gId]['ref']
                try:
                    already_included = elements[gId]
                except:
                    elements[gId] = dict()
                    elements[gId]['ref'] = genomeRef  # the key line
                    self.log(console,"adding element "+gId+" : "+genomeRef)  # DEBUG
            

        # load the method provenance from the context object
        #
        self.log(console,"SETTING PROVENANCE")  # DEBUG
        provenance = [{}]
        if 'provenance' in ctx:
            provenance = ctx['provenance']
        # add additional info to provenance here, in this case the input data object reference
        try:
            prov_defined = provenance[0]['input_ws_objects']
        except:
            provenance[0]['input_ws_objects'] = []
        provenance[0]['input_ws_objects'].append(params['workspace_name']+'/'+params['input_genome_names'])
        provenance[0]['input_ws_objects'].append(params['workspace_name']+'/'+params['input_genomeset_name'])
        provenance[0]['service'] = 'kb_util_dylan'
        provenance[0]['method'] = 'KButil_Add_Genomes_to_GenomeSet'


        # Store output object
        #
        if len(invalid_msgs) == 0:
            self.log(console,"SAVING GENOMESET")  # DEBUG
            output_GenomeSet = {
                              'description': params['desc'],
                              'elements': elements
                            }

            new_obj_info = ws.save_objects({
                            'workspace': params['workspace_name'],
                            'objects':[{
                                    'type': 'KBaseSearch.GenomeSet',
                                    'data': output_GenomeSet,
                                    'name': params['output_name'],
                                    'meta': {},
                                    'provenance': provenance
                                }]
                        })


        # build output report object
        #
        self.log(console,"BUILDING REPORT")  # DEBUG
        if len(invalid_msgs) == 0:
            self.log(console,"genomes in output set "+params['output_name']+": "+str(len(elements.keys())))
            report += 'genomes in output set '+params['output_name']+': '+str(len(elements.keys()))+"\n"
            reportObj = {
                'objects_created':[{'ref':params['workspace_name']+'/'+params['output_name'], 'description':'KButil_Add_Genomes_to_GenomeSet'}],
                'text_message':report
                }
        else:
            report += "FAILURE:\n\n"+"\n".join(invalid_msgs)+"\n"
            reportObj = {
                'objects_created':[],
                'text_message':report
                }

        reportName = 'kb_util_dylan_add_genomes_to_genomeset_report_'+str(uuid.uuid4())
        ws = workspaceService(self.workspaceURL, token=ctx['token'])
        report_obj_info = ws.save_objects({
#                'id':info[6],
                'workspace':params['workspace_name'],
                'objects':[
                    {
                        'type':'KBaseReport.Report',
                        'data':reportObj,
                        'name':reportName,
                        'meta':{},
                        'hidden':1,
                        'provenance':provenance
                    }
                ]
            })[0]


        # Build report and return
        #
        self.log(console,"BUILDING RETURN OBJECT")
        returnVal = { 'report_name': reportName,
                      'report_ref': str(report_obj_info[6]) + '/' + str(report_obj_info[0]) + '/' + str(report_obj_info[4]),
                      }
        self.log(console,"KButil_Add_Genomes_to_GenomeSet DONE")
        #END KButil_Add_Genomes_to_GenomeSet

        # At some point might do deeper type checking...
        if not isinstance(returnVal, dict):
            raise ValueError('Method KButil_Add_Genomes_to_GenomeSet return value ' +
                             'returnVal is not type dict as required.')
        # return the results
        return [returnVal]

    def KButil_Concat_MSAs(self, ctx, params):
        """
        :param params: instance of type "KButil_Concat_MSAs_Params"
           (KButil_Concat_MSAs() ** **  Method for Concatenating MSAs into a
           combined MSA) -> structure: parameter "workspace_name" of type
           "workspace_name" (** The workspace object refs are of form: ** ** 
           objects = ws.get_objects([{'ref':
           params['workspace_id']+'/'+params['obj_name']}]) ** ** "ref" means
           the entire name combining the workspace id and the object name **
           "id" is a numerical identifier of the workspace or object, and
           should just be used for workspace ** "name" is a string identifier
           of a workspace or object.  This is received from Narrative.),
           parameter "input_names" of type "data_obj_name", parameter
           "output_name" of type "data_obj_name", parameter "desc" of String,
           parameter "blanks_flag" of Long
        :returns: instance of type "KButil_Concat_MSAs_Output" -> structure:
           parameter "report_name" of type "data_obj_name", parameter
           "report_ref" of type "data_obj_ref"
        """
        # ctx is the context object
        # return variables are: returnVal
        #BEGIN KButil_Concat_MSAs
        console = []
        invalid_msgs = []
        self.log(console,'Running KButil_Concat_MSAs with params=')
        self.log(console, "\n"+pformat(params))
        report = ''
#        report = 'Running KButil_Concat_MSAs with params='
#        report += "\n"+pformat(params)


        #### do some basic checks
        #
        if 'workspace_name' not in params:
            raise ValueError('workspace_name parameter is required')
        if 'desc' not in params:
            raise ValueError('desc parameter is required')
        if 'input_names' not in params:
            raise ValueError('input_names parameter is required')
        if 'output_name' not in params:
            raise ValueError('output_name parameter is required')

        if len(params['input_names']) < 2:
            self.log(console,"Must provide more than one MSA")
            self.log(invalid_msgs,"Must provide more than one MSA")


        # Build FeatureSet
        #
        row_order = []
        alignment = {}
        curr_pos = 0
        MSA_seen = {}
        discard_set = {}
        sequence_type = None
        for MSA_i,MSA_name in enumerate(params['input_names']):
            if len(params['input_names']) < 2:  # too lazy to reindent the block
                continue

            if not MSA_name in MSA_seen.keys():
                MSA_seen[MSA_name] = True
            else:
                self.log(console,"repeat MSA_name: '"+MSA_name+"'")
                self.log(invalid_msgs,"repeat MSA_name: '"+MSA_name+"'")
                continue

            try:
                ws = workspaceService(self.workspaceURL, token=ctx['token'])
                objects = ws.get_objects([{'ref': params['workspace_name']+'/'+MSA_name}])
                data = objects[0]['data']
                info = objects[0]['info']
                # Object Info Contents
                # absolute ref = info[6] + '/' + info[0] + '/' + info[4]
                # 0 - obj_id objid
                # 1 - obj_name name
                # 2 - type_string type
                # 3 - timestamp save_date
                # 4 - int version
                # 5 - username saved_by
                # 6 - ws_id wsid
                # 7 - ws_name workspace
                # 8 - string chsum
                # 9 - int size 
                # 10 - usermeta meta
                type_name = info[2].split('.')[1].split('-')[0]

            except Exception as e:
                raise ValueError('Unable to fetch input_name object from workspace: ' + str(e))
                #to get the full stack trace: traceback.format_exc()

            if type_name != 'MSA':
                raise ValueError("Bad Type:  Should be MSA instead of '"+type_name+"'")

            this_MSA = data
            this_genomes_seen = {}

            # set sequence_type
            try:
                this_sequence_type = this_MSA['sequence_type']
                if sequence_type == None:
                    sequence_type = this_sequence_type
                elif this_sequence_type != sequence_type:
                    self.log(invalid_msgs,"inconsistent sequence type for MSA "+MSA_name+" '"+this_sequence_type+"' doesn't match '"+sequence_type+"'")
                    continue
            except:
                pass

            # build row_order
            this_row_order = []
            if 'row_order' in this_MSA.keys():
                #self.log(console,"IN row_order A")  # DEBUG
                this_row_order = this_MSA['row_order']
            else:
                #self.log(console,"IN row_order B")  # DEBUG
                this_row_order = sorted(this_MSA['alignment'].keys())

            # DEBUG
            #for row_id in this_row_order:
            #    self.log(console,"ROW_ORDER_ID: '"+row_id+"'")
            #for row_id in sorted(this_MSA['alignment']):
            #    self.log(console,"ALIGNMENT_ID: '"+row_id+"'")


            # concat alignments using base genome_id to unify (input rows are probably gene ids)
            this_aln_len = len(this_MSA['alignment'][this_row_order[0]])
            genome_row_ids_updated = {}
            for row_id in this_row_order:
                id_pieces = re.split('\.', row_id)
                if len(id_pieces) >= 2:
                    genome_id = ".".join(id_pieces[0:2])  # just want genome_id
                else:
                    genome_id = row_id

                # can't have repeat genome_ids (i.e. no paralogs allowed)
                try:
                    genome_id_seen = this_genomes_seen[genome_id]
                    self.log(console,"only one feature per genome permitted in a given MSA.  MSA: "+MSA_name+" genome_id: "+genome_id+" row_id: "+row_id)
                    self.log(invalid_msgs,"only one feature per genome permitted in a given MSA.  MSA: "+MSA_name+" genome_id: "+genome_id+" row_id: "+row_id)
                    continue
                except:
                    this_genomes_seen[genome_id] = True

                this_row_len = len(this_MSA['alignment'][row_id])
                if this_row_len != this_aln_len:
                    self.log(invalid_msgs,"inconsistent alignment len in "+MSA_name+": first_row_len="+str(this_aln_len)+" != "+str(this_row_len)+" ("+row_id+")")
                    continue

                # create new rows
                if genome_id not in alignment.keys():
                    row_order.append(genome_id)
                    alignment[genome_id] = ''
                    if MSA_i > 0:
                        #self.log(console,"NOT IN OLD MSA: "+genome_id)  # DEBUG
                        discard_set[genome_id] = True
                        alignment[genome_id] += ''.join(['-' for s in range(curr_pos)])
                #else:  # DEBUG
                    #self.log(console,"SEEN IN MSA: "+genome_id)  # DEBUG

                # add new 
                genome_row_ids_updated[genome_id] = True
                alignment[genome_id] += this_MSA['alignment'][row_id]

            # append blanks for rows that weren't in new MSA
            if MSA_i > 0:
                for genome_id in alignment.keys():
                    try:
                        updated = genome_row_ids_updated[genome_id]
                    except:
                        #self.log(console,"NOT IN NEW MSA: "+genome_id)  # DEBUG
                        discard_set[genome_id] = True
                        alignment[genome_id] += ''.join(['-' for s in range(this_aln_len)])
            # update curr_pos
            curr_pos += this_aln_len
            
            # report
            if len(invalid_msgs) == 0:
                report += 'num rows in input set '+MSA_name+': '+str(len(this_row_order))+" "+str(this_row_order)+"\n"
                self.log(console,'num rows in input set '+MSA_name+': '+str(len(this_row_order)))
                self.log(console,'row_ids in input set '+MSA_name+': '+str(this_row_order))

        # report which are incomplete rows (regardless of whether discarding)
        if len(invalid_msgs) == 0:
            for genome_id in row_order:
                try:
                    discard = discard_set[genome_id]
                    self.log(console,'incomplete row: '+genome_id+"\n")
                    report += 'incomplete row: '+genome_id
                except:
                    self.log(console,'complete row: '+genome_id+"\n")
                    report += 'complete row: '+genome_id
        

            # remove incomplete rows if not adding blanks
            if 'blanks_flag' in params and params['blanks_flag'] != None and params['blanks_flag'] == 0:
                new_row_order = []
                new_alignment = {}
                for genome_id in row_order:
                    try:
                        discard = discard_set[genome_id]
                        self.log(console,'discarding row: '+genome_id+"\n")
                        report += 'discarding row: '+genome_id
                    except:
                        new_row_order.append(genome_id)
                        new_alignment[genome_id] = alignment[genome_id]
                    row_order = new_row_order
                    alignment = new_alignment

            # report which rows are retained
            for genome_id in row_order:
                self.log(console,'output MSA contains row: '+genome_id+"\n")
                report += 'output MSA contains row: '+genome_id


        # load the method provenance from the context object
        #
        self.log(console,"SETTING PROVENANCE")  # DEBUG
        provenance = [{}]
        if 'provenance' in ctx:
            provenance = ctx['provenance']
        # add additional info to provenance here, in this case the input data object reference
        provenance[0]['input_ws_objects'] = []
        for MSA_name in params['input_names']:
            provenance[0]['input_ws_objects'].append(params['workspace_name']+'/'+MSA_name)
        provenance[0]['service'] = 'kb_util_dylan'
        provenance[0]['method'] = 'KButil_Concat_MSAs'


        # DEBUG: check alignment and row_order
        #for genome_id in row_order:
        #    self.log(console,"AFTER ROW_ORDER: "+genome_id)
        #for genome_id in alignment.keys():
        #    self.log(console,"AFTER ALIGNMENT: "+genome_id+",\t"+alignment[genome_id])


        # Store output object
        #
        if len(invalid_msgs) == 0:
            self.log(console,"SAVING OUTPUT MSA")  # DEBUG
            output_MSA = {
                       'name': params['output_name'],
                       'description': params['desc'],
                       'row_order': row_order,
                       'alignment': alignment,
                       'alignment_length': len(alignment[row_order[0]])
                     }
            if sequence_type != None:
                output_MSA['sequence_type'] = sequence_type

                new_obj_info = ws.save_objects({
                            'workspace': params['workspace_name'],
                            'objects':[{
                                    'type': 'KBaseTrees.MSA',
                                    'data': output_MSA,
                                    'name': params['output_name'],
                                    'meta': {},
                                    'provenance': provenance
                                }]
                        })


        # build output report object
        #
        self.log(console,"BUILDING REPORT")  # DEBUG
        if len(invalid_msgs) == 0:
            self.log(console,"rows in output MSA "+params['output_name']+": "+str(len(row_order)))
            report += 'rows in output MSA '+params['output_name']+': '+str(len(row_order))+"\n"
            reportObj = {
                'objects_created':[{'ref':params['workspace_name']+'/'+params['output_name'], 'description':'KButil_Concat_MSAs'}],
                'text_message':report
                }
        else:
            report += "FAILURE:\n\n"+"\n".join(invalid_msgs)+"\n"
            reportObj = {
                'objects_created':[],
                'text_message':report
                }

        reportName = 'kb_util_dylan_concat_msas_report_'+str(uuid.uuid4())
        ws = workspaceService(self.workspaceURL, token=ctx['token'])
        report_obj_info = ws.save_objects({
#                'id':info[6],
                'workspace':params['workspace_name'],
                'objects':[
                    {
                        'type':'KBaseReport.Report',
                        'data':reportObj,
                        'name':reportName,
                        'meta':{},
                        'hidden':1,
                        'provenance':provenance
                    }
                ]
            })[0]


        # Build report and return
        #
        self.log(console,"BUILDING RETURN OBJECT")
        returnVal = { 'report_name': reportName,
                      'report_ref': str(report_obj_info[6]) + '/' + str(report_obj_info[0]) + '/' + str(report_obj_info[4]),
                      }
        self.log(console,"KButil_Concat_MSAs DONE")
        #END KButil_Concat_MSAs

        # At some point might do deeper type checking...
        if not isinstance(returnVal, dict):
            raise ValueError('Method KButil_Concat_MSAs return value ' +
                             'returnVal is not type dict as required.')
        # return the results
        return [returnVal]

    def KButil_Split_Reads(self, ctx, input_params):
        """
        :param params: instance of type "KButil_Split_Reads_Params"
           (KButil_Split_Reads() ** **  Method for spliting a ReadsLibrary
           into evenly sized ReadsLibraries) -> structure: parameter
           "workspace_name" of type "workspace_name" (** The workspace object
           refs are of form: ** **    objects = ws.get_objects([{'ref':
           params['workspace_id']+'/'+params['obj_name']}]) ** ** "ref" means
           the entire name combining the workspace id and the object name **
           "id" is a numerical identifier of the workspace or object, and
           should just be used for workspace ** "name" is a string identifier
           of a workspace or object.  This is received from Narrative.),
           parameter "input_name" of type "data_obj_name", parameter
           "output_name" of type "data_obj_name", parameter "split_num" of
           Long, parameter "desc" of String
        :returns: instance of type "KButil_Split_Reads_Output" -> structure:
           parameter "report_name" of type "data_obj_name", parameter
           "report_ref" of type "data_obj_ref"
        """
        # ctx is the context object
        # return variables are: returnVal
        #BEGIN KButil_Split_Reads
        console = []
        self.log(console, 'Running KButil_Split_Reads() with parameters: ')
        self.log(console, "\n"+pformat(input_params))
        
        token = ctx['token']
        wsClient = workspaceService(self.workspaceURL, token=token)
        headers = {'Authorization': 'OAuth '+token}
        env = os.environ.copy()
        env['KB_AUTH_TOKEN'] = token
        
        SERVICE_VER = 'dev'  # DEBUG

        # param checks
        required_params = ['input_name', 
                           'output_name'
                           ]
        for required_param in required_params:
            if required_param not in input_params or input_params[required_param] == None:
                raise ValueError ("Must define required param: '"+required_param+"'")
            
        # and param defaults
        defaults = { 'split_num': 10
                   }
        for arg in defaults.keys():
            if arg not in input_params or input_params[arg] == None or input_params[arg] == '':
                input_params[arg] = defaults[arg]

        # load provenance
        provenance = [{}]
        if 'provenance' in ctx:
            provenance = ctx['provenance']
        provenance[0]['input_ws_objects']=[str(input_params['workspace_name'])+'/'+str(input_params['input_name'])]


        # Determine whether read library is of correct type
        #
        try:
            # object_info tuple
            [OBJID_I, NAME_I, TYPE_I, SAVE_DATE_I, VERSION_I, SAVED_BY_I, WSID_I, WORKSPACE_I, CHSUM_I, SIZE_I, META_I] = range(11)
            
            input_reads_ref = str(input_params['workspace_name'])+'/'+str(input_params['input_name'])
            input_reads_obj_info = wsClient.get_object_info_new ({'objects':[{'ref':input_reads_ref}]})[0]
            input_reads_obj_type = input_reads_obj_info[TYPE_I]
            #input_reads_obj_version = input_reads_obj_info[VERSION_I]  # this is object version, not type version

        except Exception as e:
            raise ValueError('Unable to get read library object from workspace: (' + str(input_params['input_reads_ref']) +')' + str(e))

        input_reads_obj_type = re.sub ('-[0-9]+\.[0-9]+$', "", input_reads_obj_type)  # remove trailing version

        acceptable_types = ["KBaseFile.PairedEndLibrary", "KBaseFile.SingleEndLibrary"]
        if input_reads_obj_type not in acceptable_types:
            raise ValueError ("Input reads of type: '"+input_reads_obj_type+"'.  Must be one of "+", ".join(acceptable_types))


        # Download Reads
        #
        try:
            readsUtils_Client = ReadsUtils (url=self.callbackURL, token=ctx['token'])  # SDK local
            
            readsLibrary = readsUtils_Client.download_reads ({'read_libraries': [input_reads_ref],
                                                             'interleaved': 'false'
                                                             })
        except Exception as e:
            raise ValueError('Unable to get read library object from workspace: (' + str(input_reads_ref) +")\n" + str(e))


        # Paired End
        #
        if input_reads_obj_type == "KBaseFile.PairedEndLibrary":

            # Download reads Libs to FASTQ files
            input_fwd_file_path = readsLibrary['files'][input_reads_ref]['files']['fwd']
            input_rev_file_path = readsLibrary['files'][input_reads_ref]['files']['rev']
            sequencing_tech     = readsLibrary['files'][input_reads_ref]['sequencing_tech']

            # file paths
            input_fwd_path = re.sub ("\.fastq$", "", input_fwd_file_path)
            input_fwd_path = re.sub ("\.FASTQ$", "", input_fwd_path)
            input_rev_path = re.sub ("\.fastq$", "", input_rev_file_path)
            input_rev_path = re.sub ("\.FASTQ$", "", input_rev_path)
            output_fwd_paired_file_path_base   = input_fwd_path+"_fwd_paired"
            output_rev_paired_file_path_base   = input_rev_path+"_rev_paired"
            output_fwd_unpaired_file_path = input_fwd_path+"_fwd_unpaired.fastq"
            output_rev_unpaired_file_path = input_rev_path+"_rev_unpaired.fastq"
            # set up for file io
            total_paired_reads = 0
            total_paired_reads_by_set = []
            total_unpaired_fwd_reads = 0
            total_unpaired_rev_reads = 0
            fwd_ids = dict()
            paired_lib_i = dict()
            unpaired_buf_size = 0
            paired_buf_size = 1000000

            # read fwd file to get fwd ids
            with open (input_fwd_file_path, 'r', 0) as input_reads_file_handle:
                for line in input_reads_file_handle:
                    if line.startswith('@'):
                        read_id = re.sub ("[ \t]+.*", "", line)
                        fwd_ids[read_id] = True

            # determine paired and unpaired rev, split paired rev, write unpaired rev, and store lib_i for paired
            paired_output_reads_file_handles = []
            for lib_i in range(input_params['split_num']):
                paired_output_reads_file_handles[lib_i] = open (output_rev_paired_file_path_base+"-"+str(lib_i)+".fastq", 'w', paired_buf_size)
                total_paired_reads_by_set[lib_i] = 0

            rec_buf = []
            unpaired_rev_buf = []
            last_read_id = None
            paired_cnt = 0
            capture_type_paired = False

            with open (input_rev_file_path, 'r', 0) as input_reads_file_handle:
                for line in input_reads_file_handle:
                    if line.startswith('@'):
                        if last_read_id != None:
                            if capture_type_paired:
                                lib_i = paired_cnt % input_params['split_num']
                                total_paired_reads_by_set[lib_i] += 1
                                paired_lib_i[last_read_id] = lib_i
                                paired_output_reads_file_handles[lib_i].writelines(rec_buf)
                                paired_cnt += 1
                            else:
                                unpaired_rev_buf.extend(rec_buf)
                            rec_buf = []
                        last_read_id = read_id = re.sub ("[ \t]+.*", "", line)
                        try:
                            found = fwd_ids[read_id]
                            total_paired_reads += 1
                            capture_type_paired = True
                        except:
                            total_unpaired_rev_reads += 1
                            capture_type_paired = False
                    rec_buf.append(line)

            for output_handle in paired_output_reads_file_handles:
                output_handle.close()

            output_reads_file_handle = open (output_rev_unpaired_file_path, 'w', 0)
            output_reads_file_handle.writelines(unpaired_rev_buf)
            output_reads_file_handle.close()


            # split fwd paired and write unpaired fwd
            paired_output_reads_file_handles = []
            for lib_i in range(input_params['split_num']):
                paired_output_reads_file_handles[lib_i] = open (output_fwd_paired_file_path_base+"-"+str(lib_i)+".fastq", 'w', paired_buf_size)

            rec_buf = []
            unpaired_fwd_buf = []
            last_read_id = None
            paired_cnt = 0
            capture_type_paired = False

            with open (input_fwd_file_path, 'r', 0) as input_reads_file_handle:
                for line in input_reads_file_handle:
                    if line.startswith('@'):
                        if last_read_id != None:
                            if capture_type_paired:
                                lib_i = paired_lib_i[last_read_id]
                                paired_output_reads_file_handles[lib_i].writelines(rec_buf)
                                paired_cnt += 1
                            else:
                                unpaired_fwd_buf.extend(rec_buf)
                            rec_buf = []
                        last_read_id = read_id = re.sub ("[ \t]+.*", "", line)
                        try:
                            found = paired_lib_i[read_id]
                            capture_type_paired = True
                        except:
                            total_unpaired_fwd_reads += 1
                            capture_type_paired = False
                    rec_buf.append(line)

            for output_handle in paired_output_reads_file_handles:
                output_handle.close()

            output_reads_file_handle = open (output_fwd_unpaired_file_path, 'w', 0)
            output_reads_file_handle.writelines(unpaired_fwd_buf)
            output_reads_file_handle.close()


            # store report
            #
            report += "TOTAL PAIRED READS: "+total_paired_reads+"\n"
            report += "TOTAL UNPAIRED FWD READS: " +total_unpaired_fwd_reads+"\n"
            report += "TOTAL UNPAIRED REV READS: " +total_unpaired_rev_reads+"\n"
            report += "\n"
            for lib_i in range(input_params['split_num']):
                report += "PAIRED READS IN SET "+str(lib_i)+": "+str(total_paired_reads_by_set[lib_i])+"\n"


            # upload paired reads
            #
            paired_obj_refs = []
            for lib_i in range(input_params['split_num']):
                output_fwd_paired_file_path = output_fwd_paired_file_path_base+"-"+str(lib_i)+".fastq"
                output_rev_paired_file_path = output_rev_paired_file_path_base+"-"+str(lib_i)+".fastq"
                if not os.path.isfile (output_fwd_paired_file_path) \
                     or os.path.getsize (output_fwd_paired_file_path) == 0 \
                   or not os.path.isfile (output_rev_paired_file_path) \
                     or os.path.getsize (output_rev_paired_file_path) == 0:
                    
                    raise ValueError ("failed to create paired output")
                else:
                    output_obj_name = input_params['output_name']+'_paired'+str(lib_i)
                    self.log(console, 'Uploading paired reads: '+output_obj_name)
                    paired_obj_refs[lib_i] = readsUtils_Client.upload_reads ({ 'wsname': str(input_params['workspace_name']),
                                                                                      'name': output_obj_name,
                                                                                      'sequencing_tech': sequencing_tech,
                                                                                      'fwd_file': output_fwd_paired_file_path,
                                                                                      'rev_file': output_rev_paired_file_path
                                                                                      })['obj_ref']
                    

            # upload reads forward unpaired
            unpaired_fwd_ref = None
            if os.path.isfile (output_fwd_unpaired_file_path) \
                and os.path.getsize (output_fwd_unpaired_file_path) != 0:

                output_obj_name = input_params['output_reads_name']+'_unpaired_fwd'
                self.log(console, '\nUploading trimmed unpaired forward reads: '+output_obj_name)
                unpaired_fwd_ref = readsUtils_Client.upload_reads ({ 'wsname': str(input_params['workspace_name']),
                                                                     'name': output_obj_name,
                                                                     'sequencing_tech': sequencing_tech,
                                                                     'fwd_file': output_fwd_unpaired_file_path
                                                                     })['obj_ref']
                

            # upload reads reverse unpaired
            unpaired_rev_ref = None
            if os.path.isfile (output_rev_unpaired_file_path) \
                and os.path.getsize (output_rev_unpaired_file_path) != 0:

                output_obj_name = input_params['output_reads_name']+'_unpaired_rev'
                self.log(console, '\nUploading trimmed unpaired reverse reads: '+output_obj_name)
                unpaired_rev_ref = readsUtils_Client.upload_reads ({ 'wsname': str(input_params['workspace_name']),
                                                                     'name': output_obj_name,
                                                                     'sequencing_tech': sequencing_tech,
                                                                     'fwd_file': output_rev_unpaired_file_path
                                                                     })['obj_ref']
                

        # SingleEndLibrary
        #
        elif input_reads_obj_type == "KBaseFile.SingleEndLibrary":
            self.log(console, "Downloading Single End reads file...")

            # Download reads Libs to FASTQ files
            input_fwd_file_path = readsLibrary['files'][input_reads_ref]['files']['fwd']
            sequencing_tech     = readsLibrary['files'][input_reads_ref]['sequencing_tech']

            # file paths
            input_fwd_path = re.sub ("\.fastq$", "", input_fwd_file_path)
            input_fwd_path = re.sub ("\.FASTQ$", "", input_fwd_path)
            output_fwd_paired_file_path_base   = input_fwd_path+"_fwd_paired"

            # set up for file io
            total_paired_reads = 0
            total_paired_reads_by_set = []
            paired_buf_size = 1000000

            # split fwd paired
            paired_output_reads_file_handles = []
            for lib_i in range(input_params['split_num']):
                paired_output_reads_file_handles[lib_i] = open (output_fwd_paired_file_path_base+"-"+str(lib_i)+".fastq", 'w', paired_buf_size)
                total_paired_reads_by_set[lib_i] = 0

            rec_buf = []
            last_read_id = None
            paired_cnt = 0

            with open (input_fwd_file_path, 'r', 0) as input_reads_file_handle:
                for line in input_reads_file_handle:
                    if line.startswith('@'):
                        if last_read_id != None:
                            lib_i = paired_cnt % input_params['split_num']
                            total_paired_reads_by_set[lib_i] += 1
                            paired_output_reads_file_handles[lib_i].writelines(rec_buf)
                            paired_cnt += 1
                            rec_buf = []
                        last_read_id = read_id = re.sub ("[ \t]+.*", "", line)
                    rec_buf.append(line)

            for output_handle in paired_output_reads_file_handles:
                output_handle.close()


            # store report
            #
            for lib_i in range(input_params['split_num']):
                report += "PAIRED READS IN SET "+str(lib_i)+": "+str(total_paired_reads_by_set[lib_i])+"\n"


            # upload paired reads
            #
            paired_obj_refs = []
            for lib_i in range(input_params['split_num']):
                output_fwd_paired_file_path = output_fwd_paired_file_path_base+"-"+str(lib_i)+".fastq"
                if not os.path.isfile (output_fwd_paired_file_path) \
                        or os.path.getsize (output_fwd_paired_file_path) == 0:
                    
                    raise ValueError ("failed to create paired output")
                else:
                    output_obj_name = input_params['output_name']+'-'+str(lib_i)
                    self.log(console, 'Uploading paired reads: '+output_obj_name)
                    paired_obj_refs[lib_i] = readsUtils_Client.upload_reads ({ 'wsname': str(input_params['workspace_name']),
                                                                               'name': output_obj_name,
                                                                               'sequencing_tech': sequencing_tech,
                                                                               'fwd_file': output_fwd_paired_file_path
                                                                               })['obj_ref']
                    
        else:
            raise ValueError ("unknown ReadLibrary type as input: "+str(input_reads_obj_type))


        # save output readsSet
        #
        items = []
        for lib_i,lib_ref in enumerate(paired_obj_refs):
            label = input_params['output_name']+'-'+str(lib_i)
            items.append({'ref': lib_ref,
                          'label': label
                          #'data_attachment': ,
                          #'info':
                              })
        description = input_params['desc']
        output_readsSet_obj = { 'description': input_params['desc'],
                                'items': items
                                }
        output_readsSet_name = str(input_params['output_name'])
        readsSet_ref = setAPI_Client.save_reads_set_v1 ({'workspace_name': input_params['workspace_name'],
                                                         'output_object_name': output_readsSet_name,
                                                         'data': output_readsSet_obj
                                                         })['set_ref']
                              

        # build report
        #
        reportObj = {'objects_created':[], 
                     'text_message': report}

        reportObj['objects_created'].append({'ref':readsSet_ref,
                                             'description':input_params['desc']})

        if unpaired_fwd_ref != None:
            reportObj['objects_created'].append({'ref':unpaired_fwd_ref,
                                                 'description':input_params['desc']+" unpaired fwd reads"})

        if unpaired_rev_ref != None:
            reportObj['objects_created'].append({'ref':unpaired_rev_ref,
                                                 'description':input_params['desc']+" unpaired rev reads"})


        # save report object
        #
        report = KBaseReport(self.callbackURL, token=ctx['token'], service_ver=SERVICE_VER)
        report_info = report.create({'report':reportObj, 'workspace_name':input_params['input_ws']})

        output = { 'report_name': report_info['name'], 'report_ref': report_info['ref'] }
        #END KButil_Split_Reads

        # At some point might do deeper type checking...
        if not isinstance(returnVal, dict):
            raise ValueError('Method KButil_Split_Reads return value ' +
                             'returnVal is not type dict as required.')
        # return the results
        return [returnVal]

    def status(self, ctx):
        #BEGIN_STATUS
        returnVal = {'state': "OK", 'message': "", 'version': self.VERSION, 
                     'git_url': self.GIT_URL, 'git_commit_hash': self.GIT_COMMIT_HASH}
        #END_STATUS
        return [returnVal]
