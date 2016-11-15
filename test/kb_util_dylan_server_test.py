import unittest
import os
import json
import time
import requests

from os import environ
from ConfigParser import ConfigParser
from requests_toolbelt import MultipartEncoder
from pprint import pprint

from biokbase.workspace.client import Workspace as workspaceService
from biokbase.AbstractHandle.Client import AbstractHandle as HandleService
from kb_util_dylan.kb_util_dylanImpl import kb_util_dylan


class kb_util_dylanTest(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        token = environ.get('KB_AUTH_TOKEN', None)
        cls.token = token
        cls.ctx = {'token': token, 'provenance': [{'service': 'kb_util_dylan',
            'method': 'please_never_use_it_in_production', 'method_params': []}],
            'authenticated': 1}
        config_file = environ.get('KB_DEPLOYMENT_CONFIG', None)
        cls.cfg = {}
        config = ConfigParser()
        config.read(config_file)
        for nameval in config.items('kb_util_dylan'):
            cls.cfg[nameval[0]] = nameval[1]
        cls.wsURL = cls.cfg['workspace-url']
        cls.wsClient = workspaceService(cls.wsURL, token=token)
        cls.serviceImpl = kb_util_dylan(cls.cfg)

        cls.shockURL = cls.cfg['shock-url']
        cls.handleURL = cls.cfg['handle-service-url']


    @classmethod
    def tearDownClass(cls):
        if hasattr(cls, 'wsName'):
            cls.wsClient.delete_workspace({'workspace': cls.wsName})
            print('Test workspace was deleted')

    def getWsClient(self):
        return self.__class__.wsClient

    def getWsName(self):
        if hasattr(self.__class__, 'wsName'):
            return self.__class__.wsName
        suffix = int(time.time() * 1000)
        wsName = "test_kb_util_dylan_" + str(suffix)
        ret = self.getWsClient().create_workspace({'workspace': wsName})
        self.__class__.wsName = wsName
        return wsName

    # call this method to get the WS object info of a Paired End Library (will
    # upload the example data if this is the first time the method is called during tests)
    def getPairedEndLibInfo(self):
        if hasattr(self.__class__, 'pairedEndLibInfo'):
            return self.__class__.pairedEndLibInfo
        # 1) upload files to shock
        token = self.ctx['token']
        forward_shock_file = self.upload_file_to_shock('data/small.forward.fq')
        reverse_shock_file = self.upload_file_to_shock('data/small.reverse.fq')
        #pprint(forward_shock_file)
        #pprint(reverse_shock_file)

        # 2) create handle
        hs = HandleService(url=self.handleURL, token=token)
        forward_handle = hs.persist_handle({
                                        'id' : forward_shock_file['id'], 
                                        'type' : 'shock',
                                        'url' : self.shockURL,
                                        'file_name': forward_shock_file['file']['name'],
                                        'remote_md5': forward_shock_file['file']['checksum']['md5']}
)

        reverse_handle = hs.persist_handle({
                                        'id' : reverse_shock_file['id'], 
                                        'type' : 'shock',
                                        'url' : self.shockURL,
                                        'file_name': reverse_shock_file['file']['name'],
                                        'remote_md5': reverse_shock_file['file']['checksum']['md5']})

        # 3) save to WS
        paired_end_library = {
            'lib1': {
                'file': {
                    'hid':forward_handle,
                    'file_name': forward_shock_file['file']['name'],
                    'id': forward_shock_file['id'],
                    'url': self.shockURL,
                    'type':'shock',
                    'remote_md5':forward_shock_file['file']['checksum']['md5']
                },
                'encoding':'UTF8',
                'type':'fastq',
                'size':forward_shock_file['file']['size']
            },
            'lib2': {
                'file': {
                    'hid':reverse_handle,
                    'file_name': reverse_shock_file['file']['name'],
                    'id': reverse_shock_file['id'],
                    'url': self.shockURL,
                    'type':'shock',
                    'remote_md5':reverse_shock_file['file']['checksum']['md5']
                },
                'encoding':'UTF8',
                'type':'fastq',
                'size':reverse_shock_file['file']['size']

            },
            'interleaved':0,
            'sequencing_tech':'artificial reads'
        }

        new_obj_info = self.wsClient.save_objects({
                        'workspace':self.getWsName(),
                        'objects':[
                            {
                                'type':'KBaseFile.PairedEndLibrary',
                                'data':paired_end_library,
                                'name':'test.pe.reads',
                                'meta':{},
                                'provenance':[
                                    {
                                        'service':'MEGAHIT',
                                        'method':'test_megahit'
                                    }
                                ]
                            }]
                        })
        self.__class__.pairedEndLibInfo = new_obj_info[0]
        return new_obj_info[0]




    def getImpl(self):
        return self.__class__.serviceImpl

    def getContext(self):
        return self.__class__.ctx


    def test_KButil_Split_Reads(self):
        # Prepare test objects in workspace if needed using 
        # self.getWsClient().save_objects({'workspace': self.getWsName(), 'objects': []})
        #
        # Run your method by
        # ret = self.getImpl().your_method(self.getContext(), parameters...)
        #
        # Check returned data with
        # self.assertEqual(ret[...], ...) or other unittest methods

        # figure out where the test data lives
        pe_lib_info = self.getPairedEndLibInfo()
        pprint(pe_lib_info)

        # Object Info Contents
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


        # run Split_Reads
        split_num = 4
        params = {
            'workspace_name': pe_lib_info[7],
            'input_name': pe_lib_info[1],
            'split_num': split_num,
            'output_name': 'output',
            'desc':'test split',
            #'k_min':31,
            #'k_max':51,
            #'k_step':10,
            #'k_list':[31,41],
            #'min_contig_length':199
        }

        result = self.getImpl().KButil_Split_Reads(self.getContext(),params)
        print('RESULT:')
        pprint(result)

        # check the output
        output_name = 'output_paired-0'
        output_type = 'KBaseFile.PairedEndLibrary'
        info_list = self.wsClient.get_object_info([{'ref':pe_lib_info[7] + '/' + output_name}], 1)
        self.assertEqual(len(info_list),1)
        readsLib_info = info_list[0]
        self.assertEqual(readsLib_info[1],output_name)
        self.assertEqual(readsLib_info[2].split('-')[0],output_type)

        output_name = 'output_paired-'+str(split_num)
        output_type = 'KBaseFile.PairedEndLibrary'
        info_list = self.wsClient.get_object_info([{'ref':pe_lib_info[7] + '/' + output_name}], 1)
        self.assertEqual(len(info_list),1)
        readsLib_info = info_list[0]
        self.assertEqual(readsLib_info[1],output_name)
        self.assertEqual(readsLib_info[2].split('-')[0],output_type)

        output_name = 'output_unpaired-fwd'
        output_type = 'KBaseFile.SingleEndLibrary'
        info_list = self.wsClient.get_object_info([{'ref':pe_lib_info[7] + '/' + output_name}], 1)
        self.assertEqual(len(info_list),1)
        readsLib_info = info_list[0]
        self.assertEqual(readsLib_info[1],output_name)
        self.assertEqual(readsLib_info[2].split('-')[0],output_type)

        output_name = 'output_unpaired-rev'
        output_type = 'KBaseFile.SingleEndLibrary'
        info_list = self.wsClient.get_object_info([{'ref':pe_lib_info[7] + '/' + output_name}], 1)
        self.assertEqual(len(info_list),1)
        readsLib_info = info_list[0]
        self.assertEqual(readsLib_info[1],output_name)
        self.assertEqual(readsLib_info[2].split('-')[0],output_type)

        
    @classmethod
    def upload_file_to_shock(cls, file_path):
        """
        Use HTTP multi-part POST to save a file to a SHOCK instance.
        """

        header = dict()
        header["Authorization"] = "Oauth {0}".format(cls.token)

        if file_path is None:
            raise Exception("No file given for upload to SHOCK!")

        with open(os.path.abspath(file_path), 'rb') as dataFile:
            files = {'upload': dataFile}
            response = requests.post(
                cls.shockURL + '/node', headers=header, files=files,
                stream=True, allow_redirects=True, timeout=30)

        if not response.ok:
            response.raise_for_status()

        result = response.json()

        if result['error']:
            raise Exception(result['error'][0])
        else:
            return result["data"]
