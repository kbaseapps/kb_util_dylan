#BEGIN_HEADER
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
    #BEGIN_CLASS_HEADER
    #END_CLASS_HEADER

    # config contains contents of config file in a hash or None if it couldn't
    # be found
    def __init__(self, config):
        #BEGIN_CONSTRUCTOR
        #END_CONSTRUCTOR
        pass

    def KButil_Insert_SingleEndLibrary(self, ctx, params):
        # ctx is the context object
        # return variables are: returnVal
        #BEGIN KButil_Insert_SingleEndLibrary
        #END KButil_Insert_SingleEndLibrary

        # At some point might do deeper type checking...
        if not isinstance(returnVal, dict):
            raise ValueError('Method KButil_Insert_SingleEndLibrary return value ' +
                             'returnVal is not type dict as required.')
        # return the results
        return [returnVal]

    def KButil_FASTQ_to_FASTA(self, ctx, params):
        # ctx is the context object
        # return variables are: returnVal
        #BEGIN KButil_FASTQ_to_FASTA
        #END KButil_FASTQ_to_FASTA

        # At some point might do deeper type checking...
        if not isinstance(returnVal, dict):
            raise ValueError('Method KButil_FASTQ_to_FASTA return value ' +
                             'returnVal is not type dict as required.')
        # return the results
        return [returnVal]
