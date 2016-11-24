package us.kbase.kbutildylan;

import com.fasterxml.jackson.core.type.TypeReference;
import java.io.File;
import java.io.IOException;
import java.net.URL;
import java.util.ArrayList;
import java.util.List;
import java.util.Map;
import us.kbase.auth.AuthToken;
import us.kbase.common.service.JsonClientCaller;
import us.kbase.common.service.JsonClientException;
import us.kbase.common.service.RpcContext;
import us.kbase.common.service.UnauthorizedException;

/**
 * <p>Original spec-file module name: kb_util_dylan</p>
 * <pre>
 * ** A KBase module: kb_util_dylan
 * **
 * ** This module contains basic utilities
 * </pre>
 */
public class KbUtilDylanClient {
    private JsonClientCaller caller;
    private String serviceVersion = null;


    /** Constructs a client with a custom URL and no user credentials.
     * @param url the URL of the service.
     */
    public KbUtilDylanClient(URL url) {
        caller = new JsonClientCaller(url);
    }
    /** Constructs a client with a custom URL.
     * @param url the URL of the service.
     * @param token the user's authorization token.
     * @throws UnauthorizedException if the token is not valid.
     * @throws IOException if an IOException occurs when checking the token's
     * validity.
     */
    public KbUtilDylanClient(URL url, AuthToken token) throws UnauthorizedException, IOException {
        caller = new JsonClientCaller(url, token);
    }

    /** Constructs a client with a custom URL.
     * @param url the URL of the service.
     * @param user the user name.
     * @param password the password for the user name.
     * @throws UnauthorizedException if the credentials are not valid.
     * @throws IOException if an IOException occurs when checking the user's
     * credentials.
     */
    public KbUtilDylanClient(URL url, String user, String password) throws UnauthorizedException, IOException {
        caller = new JsonClientCaller(url, user, password);
    }

    /** Constructs a client with a custom URL
     * and a custom authorization service URL.
     * @param url the URL of the service.
     * @param user the user name.
     * @param password the password for the user name.
     * @param auth the URL of the authorization server.
     * @throws UnauthorizedException if the credentials are not valid.
     * @throws IOException if an IOException occurs when checking the user's
     * credentials.
     */
    public KbUtilDylanClient(URL url, String user, String password, URL auth) throws UnauthorizedException, IOException {
        caller = new JsonClientCaller(url, user, password, auth);
    }

    /** Get the token this client uses to communicate with the server.
     * @return the authorization token.
     */
    public AuthToken getToken() {
        return caller.getToken();
    }

    /** Get the URL of the service with which this client communicates.
     * @return the service URL.
     */
    public URL getURL() {
        return caller.getURL();
    }

    /** Set the timeout between establishing a connection to a server and
     * receiving a response. A value of zero or null implies no timeout.
     * @param milliseconds the milliseconds to wait before timing out when
     * attempting to read from a server.
     */
    public void setConnectionReadTimeOut(Integer milliseconds) {
        this.caller.setConnectionReadTimeOut(milliseconds);
    }

    /** Check if this client allows insecure http (vs https) connections.
     * @return true if insecure connections are allowed.
     */
    public boolean isInsecureHttpConnectionAllowed() {
        return caller.isInsecureHttpConnectionAllowed();
    }

    /** Deprecated. Use isInsecureHttpConnectionAllowed().
     * @deprecated
     */
    public boolean isAuthAllowedForHttp() {
        return caller.isAuthAllowedForHttp();
    }

    /** Set whether insecure http (vs https) connections should be allowed by
     * this client.
     * @param allowed true to allow insecure connections. Default false
     */
    public void setIsInsecureHttpConnectionAllowed(boolean allowed) {
        caller.setInsecureHttpConnectionAllowed(allowed);
    }

    /** Deprecated. Use setIsInsecureHttpConnectionAllowed().
     * @deprecated
     */
    public void setAuthAllowedForHttp(boolean isAuthAllowedForHttp) {
        caller.setAuthAllowedForHttp(isAuthAllowedForHttp);
    }

    /** Set whether all SSL certificates, including self-signed certificates,
     * should be trusted.
     * @param trustAll true to trust all certificates. Default false.
     */
    public void setAllSSLCertificatesTrusted(final boolean trustAll) {
        caller.setAllSSLCertificatesTrusted(trustAll);
    }
    
    /** Check if this client trusts all SSL certificates, including
     * self-signed certificates.
     * @return true if all certificates are trusted.
     */
    public boolean isAllSSLCertificatesTrusted() {
        return caller.isAllSSLCertificatesTrusted();
    }
    /** Sets streaming mode on. In this case, the data will be streamed to
     * the server in chunks as it is read from disk rather than buffered in
     * memory. Many servers are not compatible with this feature.
     * @param streamRequest true to set streaming mode on, false otherwise.
     */
    public void setStreamingModeOn(boolean streamRequest) {
        caller.setStreamingModeOn(streamRequest);
    }

    /** Returns true if streaming mode is on.
     * @return true if streaming mode is on.
     */
    public boolean isStreamingModeOn() {
        return caller.isStreamingModeOn();
    }

    public void _setFileForNextRpcResponse(File f) {
        caller.setFileForNextRpcResponse(f);
    }

    public String getServiceVersion() {
        return this.serviceVersion;
    }

    public void setServiceVersion(String newValue) {
        this.serviceVersion = newValue;
    }

    /**
     * <p>Original spec-file function name: KButil_FASTQ_to_FASTA</p>
     * <pre>
     * </pre>
     * @param   params   instance of type {@link us.kbase.kbutildylan.KButilFASTQToFASTAParams KButilFASTQToFASTAParams} (original type "KButil_FASTQ_to_FASTA_Params")
     * @return   instance of type {@link us.kbase.kbutildylan.KButilFASTQToFASTAOutput KButilFASTQToFASTAOutput} (original type "KButil_FASTQ_to_FASTA_Output")
     * @throws IOException if an IO exception occurs
     * @throws JsonClientException if a JSON RPC exception occurs
     */
    public KButilFASTQToFASTAOutput kButilFASTQToFASTA(KButilFASTQToFASTAParams params, RpcContext... jsonRpcContext) throws IOException, JsonClientException {
        List<Object> args = new ArrayList<Object>();
        args.add(params);
        TypeReference<List<KButilFASTQToFASTAOutput>> retType = new TypeReference<List<KButilFASTQToFASTAOutput>>() {};
        List<KButilFASTQToFASTAOutput> res = caller.jsonrpcCall("kb_util_dylan.KButil_FASTQ_to_FASTA", args, retType, true, true, jsonRpcContext, this.serviceVersion);
        return res.get(0);
    }

    /**
     * <p>Original spec-file function name: KButil_Merge_FeatureSet_Collection</p>
     * <pre>
     * </pre>
     * @param   params   instance of type {@link us.kbase.kbutildylan.KButilMergeFeatureSetCollectionParams KButilMergeFeatureSetCollectionParams} (original type "KButil_Merge_FeatureSet_Collection_Params")
     * @return   instance of type {@link us.kbase.kbutildylan.KButilMergeFeatureSetCollectionOutput KButilMergeFeatureSetCollectionOutput} (original type "KButil_Merge_FeatureSet_Collection_Output")
     * @throws IOException if an IO exception occurs
     * @throws JsonClientException if a JSON RPC exception occurs
     */
    public KButilMergeFeatureSetCollectionOutput kButilMergeFeatureSetCollection(KButilMergeFeatureSetCollectionParams params, RpcContext... jsonRpcContext) throws IOException, JsonClientException {
        List<Object> args = new ArrayList<Object>();
        args.add(params);
        TypeReference<List<KButilMergeFeatureSetCollectionOutput>> retType = new TypeReference<List<KButilMergeFeatureSetCollectionOutput>>() {};
        List<KButilMergeFeatureSetCollectionOutput> res = caller.jsonrpcCall("kb_util_dylan.KButil_Merge_FeatureSet_Collection", args, retType, true, true, jsonRpcContext, this.serviceVersion);
        return res.get(0);
    }

    /**
     * <p>Original spec-file function name: KButil_Merge_GenomeSets</p>
     * <pre>
     * </pre>
     * @param   params   instance of type {@link us.kbase.kbutildylan.KButilMergeGenomeSetsParams KButilMergeGenomeSetsParams} (original type "KButil_Merge_GenomeSets_Params")
     * @return   instance of type {@link us.kbase.kbutildylan.KButilMergeGenomeSetsOutput KButilMergeGenomeSetsOutput} (original type "KButil_Merge_GenomeSets_Output")
     * @throws IOException if an IO exception occurs
     * @throws JsonClientException if a JSON RPC exception occurs
     */
    public KButilMergeGenomeSetsOutput kButilMergeGenomeSets(KButilMergeGenomeSetsParams params, RpcContext... jsonRpcContext) throws IOException, JsonClientException {
        List<Object> args = new ArrayList<Object>();
        args.add(params);
        TypeReference<List<KButilMergeGenomeSetsOutput>> retType = new TypeReference<List<KButilMergeGenomeSetsOutput>>() {};
        List<KButilMergeGenomeSetsOutput> res = caller.jsonrpcCall("kb_util_dylan.KButil_Merge_GenomeSets", args, retType, true, true, jsonRpcContext, this.serviceVersion);
        return res.get(0);
    }

    /**
     * <p>Original spec-file function name: KButil_Build_GenomeSet</p>
     * <pre>
     * </pre>
     * @param   params   instance of type {@link us.kbase.kbutildylan.KButilBuildGenomeSetParams KButilBuildGenomeSetParams} (original type "KButil_Build_GenomeSet_Params")
     * @return   instance of type {@link us.kbase.kbutildylan.KButilBuildGenomeSetOutput KButilBuildGenomeSetOutput} (original type "KButil_Build_GenomeSet_Output")
     * @throws IOException if an IO exception occurs
     * @throws JsonClientException if a JSON RPC exception occurs
     */
    public KButilBuildGenomeSetOutput kButilBuildGenomeSet(KButilBuildGenomeSetParams params, RpcContext... jsonRpcContext) throws IOException, JsonClientException {
        List<Object> args = new ArrayList<Object>();
        args.add(params);
        TypeReference<List<KButilBuildGenomeSetOutput>> retType = new TypeReference<List<KButilBuildGenomeSetOutput>>() {};
        List<KButilBuildGenomeSetOutput> res = caller.jsonrpcCall("kb_util_dylan.KButil_Build_GenomeSet", args, retType, true, true, jsonRpcContext, this.serviceVersion);
        return res.get(0);
    }

    /**
     * <p>Original spec-file function name: KButil_Build_GenomeSet_from_FeatureSet</p>
     * <pre>
     * </pre>
     * @param   params   instance of type {@link us.kbase.kbutildylan.KButilBuildGenomeSetFromFeatureSetParams KButilBuildGenomeSetFromFeatureSetParams} (original type "KButil_Build_GenomeSet_from_FeatureSet_Params")
     * @return   instance of type {@link us.kbase.kbutildylan.KButilBuildGenomeSetFromFeatureSetOutput KButilBuildGenomeSetFromFeatureSetOutput} (original type "KButil_Build_GenomeSet_from_FeatureSet_Output")
     * @throws IOException if an IO exception occurs
     * @throws JsonClientException if a JSON RPC exception occurs
     */
    public KButilBuildGenomeSetFromFeatureSetOutput kButilBuildGenomeSetFromFeatureSet(KButilBuildGenomeSetFromFeatureSetParams params, RpcContext... jsonRpcContext) throws IOException, JsonClientException {
        List<Object> args = new ArrayList<Object>();
        args.add(params);
        TypeReference<List<KButilBuildGenomeSetFromFeatureSetOutput>> retType = new TypeReference<List<KButilBuildGenomeSetFromFeatureSetOutput>>() {};
        List<KButilBuildGenomeSetFromFeatureSetOutput> res = caller.jsonrpcCall("kb_util_dylan.KButil_Build_GenomeSet_from_FeatureSet", args, retType, true, true, jsonRpcContext, this.serviceVersion);
        return res.get(0);
    }

    /**
     * <p>Original spec-file function name: KButil_Add_Genomes_to_GenomeSet</p>
     * <pre>
     * </pre>
     * @param   params   instance of type {@link us.kbase.kbutildylan.KButilAddGenomesToGenomeSetParams KButilAddGenomesToGenomeSetParams} (original type "KButil_Add_Genomes_to_GenomeSet_Params")
     * @return   instance of type {@link us.kbase.kbutildylan.KButilAddGenomesToGenomeSetOutput KButilAddGenomesToGenomeSetOutput} (original type "KButil_Add_Genomes_to_GenomeSet_Output")
     * @throws IOException if an IO exception occurs
     * @throws JsonClientException if a JSON RPC exception occurs
     */
    public KButilAddGenomesToGenomeSetOutput kButilAddGenomesToGenomeSet(KButilAddGenomesToGenomeSetParams params, RpcContext... jsonRpcContext) throws IOException, JsonClientException {
        List<Object> args = new ArrayList<Object>();
        args.add(params);
        TypeReference<List<KButilAddGenomesToGenomeSetOutput>> retType = new TypeReference<List<KButilAddGenomesToGenomeSetOutput>>() {};
        List<KButilAddGenomesToGenomeSetOutput> res = caller.jsonrpcCall("kb_util_dylan.KButil_Add_Genomes_to_GenomeSet", args, retType, true, true, jsonRpcContext, this.serviceVersion);
        return res.get(0);
    }

    /**
     * <p>Original spec-file function name: KButil_Concat_MSAs</p>
     * <pre>
     * </pre>
     * @param   params   instance of type {@link us.kbase.kbutildylan.KButilConcatMSAsParams KButilConcatMSAsParams} (original type "KButil_Concat_MSAs_Params")
     * @return   instance of type {@link us.kbase.kbutildylan.KButilConcatMSAsOutput KButilConcatMSAsOutput} (original type "KButil_Concat_MSAs_Output")
     * @throws IOException if an IO exception occurs
     * @throws JsonClientException if a JSON RPC exception occurs
     */
    public KButilConcatMSAsOutput kButilConcatMSAs(KButilConcatMSAsParams params, RpcContext... jsonRpcContext) throws IOException, JsonClientException {
        List<Object> args = new ArrayList<Object>();
        args.add(params);
        TypeReference<List<KButilConcatMSAsOutput>> retType = new TypeReference<List<KButilConcatMSAsOutput>>() {};
        List<KButilConcatMSAsOutput> res = caller.jsonrpcCall("kb_util_dylan.KButil_Concat_MSAs", args, retType, true, true, jsonRpcContext, this.serviceVersion);
        return res.get(0);
    }

    /**
     * <p>Original spec-file function name: KButil_Split_Reads</p>
     * <pre>
     * </pre>
     * @param   params   instance of type {@link us.kbase.kbutildylan.KButilSplitReadsParams KButilSplitReadsParams} (original type "KButil_Split_Reads_Params")
     * @return   instance of type {@link us.kbase.kbutildylan.KButilSplitReadsOutput KButilSplitReadsOutput} (original type "KButil_Split_Reads_Output")
     * @throws IOException if an IO exception occurs
     * @throws JsonClientException if a JSON RPC exception occurs
     */
    public KButilSplitReadsOutput kButilSplitReads(KButilSplitReadsParams params, RpcContext... jsonRpcContext) throws IOException, JsonClientException {
        List<Object> args = new ArrayList<Object>();
        args.add(params);
        TypeReference<List<KButilSplitReadsOutput>> retType = new TypeReference<List<KButilSplitReadsOutput>>() {};
        List<KButilSplitReadsOutput> res = caller.jsonrpcCall("kb_util_dylan.KButil_Split_Reads", args, retType, true, true, jsonRpcContext, this.serviceVersion);
        return res.get(0);
    }

    /**
     * <p>Original spec-file function name: KButil_Random_Subsample_Reads</p>
     * <pre>
     * </pre>
     * @param   params   instance of type {@link us.kbase.kbutildylan.KButilRandomSubsampleReadsParams KButilRandomSubsampleReadsParams} (original type "KButil_Random_Subsample_Reads_Params")
     * @return   instance of type {@link us.kbase.kbutildylan.KButilRandomSubsampleReadsOutput KButilRandomSubsampleReadsOutput} (original type "KButil_Random_Subsample_Reads_Output")
     * @throws IOException if an IO exception occurs
     * @throws JsonClientException if a JSON RPC exception occurs
     */
    public KButilRandomSubsampleReadsOutput kButilRandomSubsampleReads(KButilRandomSubsampleReadsParams params, RpcContext... jsonRpcContext) throws IOException, JsonClientException {
        List<Object> args = new ArrayList<Object>();
        args.add(params);
        TypeReference<List<KButilRandomSubsampleReadsOutput>> retType = new TypeReference<List<KButilRandomSubsampleReadsOutput>>() {};
        List<KButilRandomSubsampleReadsOutput> res = caller.jsonrpcCall("kb_util_dylan.KButil_Random_Subsample_Reads", args, retType, true, true, jsonRpcContext, this.serviceVersion);
        return res.get(0);
    }

    /**
     * <p>Original spec-file function name: KButil_Merge_ReadsSet_to_OneLibrary</p>
     * <pre>
     * </pre>
     * @param   params   instance of type {@link us.kbase.kbutildylan.KButilMergeReadsSetToOneLibraryParams KButilMergeReadsSetToOneLibraryParams} (original type "KButil_Merge_ReadsSet_to_OneLibrary_Params")
     * @return   instance of type {@link us.kbase.kbutildylan.KButilMergeReadsSetToOneLibraryOutput KButilMergeReadsSetToOneLibraryOutput} (original type "KButil_Merge_ReadsSet_to_OneLibrary_Output")
     * @throws IOException if an IO exception occurs
     * @throws JsonClientException if a JSON RPC exception occurs
     */
    public KButilMergeReadsSetToOneLibraryOutput kButilMergeReadsSetToOneLibrary(KButilMergeReadsSetToOneLibraryParams params, RpcContext... jsonRpcContext) throws IOException, JsonClientException {
        List<Object> args = new ArrayList<Object>();
        args.add(params);
        TypeReference<List<KButilMergeReadsSetToOneLibraryOutput>> retType = new TypeReference<List<KButilMergeReadsSetToOneLibraryOutput>>() {};
        List<KButilMergeReadsSetToOneLibraryOutput> res = caller.jsonrpcCall("kb_util_dylan.KButil_Merge_ReadsSet_to_OneLibrary", args, retType, true, true, jsonRpcContext, this.serviceVersion);
        return res.get(0);
    }

    /**
     * <p>Original spec-file function name: KButil_Merge_MultipleReadsSets_to_OneReadsSet</p>
     * <pre>
     * </pre>
     * @param   params   instance of type {@link us.kbase.kbutildylan.KButilMergeMultipleReadsSetsToOneReadsSetParams KButilMergeMultipleReadsSetsToOneReadsSetParams} (original type "KButil_Merge_MultipleReadsSets_to_OneReadsSet_Params")
     * @return   instance of type {@link us.kbase.kbutildylan.KButilMergeMultipleReadsSetsToOneReadsSetOutput KButilMergeMultipleReadsSetsToOneReadsSetOutput} (original type "KButil_Merge_MultipleReadsSets_to_OneReadsSet_Output")
     * @throws IOException if an IO exception occurs
     * @throws JsonClientException if a JSON RPC exception occurs
     */
    public KButilMergeMultipleReadsSetsToOneReadsSetOutput kButilMergeMultipleReadsSetsToOneReadsSet(KButilMergeMultipleReadsSetsToOneReadsSetParams params, RpcContext... jsonRpcContext) throws IOException, JsonClientException {
        List<Object> args = new ArrayList<Object>();
        args.add(params);
        TypeReference<List<KButilMergeMultipleReadsSetsToOneReadsSetOutput>> retType = new TypeReference<List<KButilMergeMultipleReadsSetsToOneReadsSetOutput>>() {};
        List<KButilMergeMultipleReadsSetsToOneReadsSetOutput> res = caller.jsonrpcCall("kb_util_dylan.KButil_Merge_MultipleReadsSets_to_OneReadsSet", args, retType, true, true, jsonRpcContext, this.serviceVersion);
        return res.get(0);
    }

    /**
     * <p>Original spec-file function name: KButil_Remove_Unpaired_Reads_and_Synchronize_Pairs</p>
     * <pre>
     * </pre>
     * @param   params   instance of type {@link us.kbase.kbutildylan.KButilRemoveUnpairedReadsParams KButilRemoveUnpairedReadsParams} (original type "KButil_Remove_Unpaired_Reads_Params")
     * @return   instance of type {@link us.kbase.kbutildylan.KButilRemoveUnpairedReadsOutput KButilRemoveUnpairedReadsOutput} (original type "KButil_Remove_Unpaired_Reads_Output")
     * @throws IOException if an IO exception occurs
     * @throws JsonClientException if a JSON RPC exception occurs
     */
    public KButilRemoveUnpairedReadsOutput kButilRemoveUnpairedReadsAndSynchronizePairs(KButilRemoveUnpairedReadsParams params, RpcContext... jsonRpcContext) throws IOException, JsonClientException {
        List<Object> args = new ArrayList<Object>();
        args.add(params);
        TypeReference<List<KButilRemoveUnpairedReadsOutput>> retType = new TypeReference<List<KButilRemoveUnpairedReadsOutput>>() {};
        List<KButilRemoveUnpairedReadsOutput> res = caller.jsonrpcCall("kb_util_dylan.KButil_Remove_Unpaired_Reads_and_Synchronize_Pairs", args, retType, true, true, jsonRpcContext, this.serviceVersion);
        return res.get(0);
    }

    public Map<String, Object> status(RpcContext... jsonRpcContext) throws IOException, JsonClientException {
        List<Object> args = new ArrayList<Object>();
        TypeReference<List<Map<String, Object>>> retType = new TypeReference<List<Map<String, Object>>>() {};
        List<Map<String, Object>> res = caller.jsonrpcCall("kb_util_dylan.status", args, retType, true, false, jsonRpcContext, this.serviceVersion);
        return res.get(0);
    }
}
