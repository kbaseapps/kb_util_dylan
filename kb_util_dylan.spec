/*
** A KBase module: kb_util_dylan
**
** This module contains basic utilities
*/

module kb_util_dylan {

    /* 
    ** The workspace object refs are of form:
    **
    **    objects = ws.get_objects([{'ref': params['workspace_id']+'/'+params['obj_name']}])
    **
    ** "ref" means the entire name combining the workspace id and the object name
    ** "id" is a numerical identifier of the workspace or object, and should just be used for workspace
    ** "name" is a string identifier of a workspace or object.  This is received from Narrative.
    */
    typedef string workspace_name;
    typedef string sequence;
    typedef string data_obj_name;
    typedef string data_obj_ref;


    /* KButil_Insert_SingleEndLibrary Input Params
    */
    typedef structure {
        workspace_name workspace_name;
	sequence       input_sequence;
        data_obj_name  output_name;
    } KButil_Insert_SingleEndLibrary_Params;


    /* KButil_Insert_SingleEndLibrary Output
    */
    typedef structure {
	data_obj_name report_name;
	data_obj_ref  report_ref;
/*       data_obj_ref  output_filtered_ref;
*
*        int n_initial_seqs;
*        int n_seqs_matched;
*        int n_seqs_notmatched;
*/
    } KButil_Insert_SingleEndLibrary_Output;
	

    /* KButil_FASTQ_to_FASTA Input Params
    */
    typedef structure {
        workspace_name workspace_name;
	data_obj_name  input_name;
        data_obj_name  output_name;
    } KButil_FASTQ_to_FASTA_Params;


    /* KButil_FASTQ_to_FASTA Output
    */
    typedef structure {
	data_obj_name report_name;
	data_obj_ref  report_ref;
/*       data_obj_ref  output_filtered_ref;
*
*        int n_initial_seqs;
*        int n_seqs_matched;
*        int n_seqs_notmatched;
*/
    } KButil_FASTQ_to_FASTA_Output;
	

    /* KButil_Build_FeatureSet_Collection Input Params
    */
    typedef structure {
        workspace_name workspace_name;
	data_obj_name  input_names;
        data_obj_name  output_name;
	string         desc;
    } KButil_Build_FeatureSet_Collection_Params;


    /* KButil_Build_FeatureSet_Collection Output
    */
    typedef structure {
	data_obj_name report_name;
	data_obj_ref  report_ref;
    } KButil_Build_FeatureSet_Collection_Output;
	


    /*  Method for Inserting a textarea field with FASTA or FASTQ into a SingleEndLibrary object and importing into SHOCK and WS
    */
    funcdef KButil_Insert_SingleEndLibrary (KButil_Insert_SingleEndLibrary_Params params)  returns (KButil_Insert_SingleEndLibrary_Output) authentication required;


    /*  Method for Converting a FASTQ SingleEndLibrary to a FASTA SingleEndLibrary
    */
    funcdef KButil_FASTQ_to_FASTA (KButil_FASTQ_to_FASTA_Params params)  returns (KButil_FASTQ_to_FASTA_Output) authentication required;


    /*  Method for building individual features into a FeatureSet
    */
    funcdef KButil_Build_FeatureSet_Collection (KButil_Build_FeatureSet_Collection_Params params)  returns (KButil_Build_FeatureSet_Collection_Output) authentication required;
};
