#!/usr/bin/env nextflow

// clear && nextflow run lfq.nf -ansi-log false -resume

params.input = "$baseDir/data/*.mzML"
//params.input = "/Users/lars/Code/rnaseq-nf/data/*.mzML"
params.output = "$baseDir/results"
params.database = "$baseDir/data/W82_soybase_a2v1_and_pMOZ52.fasta"
params.msgfplus = "/Volumes/GoogleDrive/Shared drives/LCMS/Analyses/soybean analysis 20210527/MSGFPlus.jar"

/*
 * channels for peptide identification and peptide quantification
 * The channels contain pairs of an unique ID (base name of the mzML file) and the mzML file itself.
 * Note that the unique ID is passed through (nearly) all processes of the workflow in order to match corresponding files
 * at later stages. For example, mapping peptide IDs to peptide features in process 'id_mapping'.
 */
Channel.fromPath(params.input).map { file -> [file.baseName, file ]}.into { ch_1; ch_2 }

/*
 * channels for the protein database and the MS-GF+ search engine executable
 */
ch_database = Channel.value(params.database)
ch_msgfplus = Channel.value(params.msgfplus)

/*
 * detect peptide features in MS1 spectral data
 */
process peptide_feature_detection {

  tag "$uid"

  input:
    tuple uid, path(file_mzML) from ch_1

  output:
    tuple uid, path("${uid}.featureXML") into ch_3

  script:
  """
  FeatureFinderMultiplex -in ${file_mzML} \\
                         -out ${uid}.featureXML \\
                         -algorithm:labels [] \\
                         -algorithm:charge 1:5 \\
                         -algorithm:rt_typical 80.0 \\
                         -algorithm:mz_tolerance 9 \\
                         -algorithm:mz_unit ppm \\
                         -algorithm:intensity_cutoff 1000.0
  """
}

ch_3.into { ch_4; ch_5 }

/*
 * precursor correction
 * The mass spectrometry machine might have fragmented an incorrect isotopic peak.
 * In that case, we shift the precursor position back to the mono-isotopic peak.
 * This correction simplifies the subsequent database search.
 */
process precursor_correction {

  tag "$uid"

  input:
    tuple uid, path(file_mzML), path(file_featureXML) from ch_2.join(ch_4)

  output:
    tuple uid, path("${uid}_corrected.mzML") into ch_6

  script:
  """
  HighResPrecursorMassCorrector -in ${file_mzML} \\
                                -feature:in ${file_featureXML} \\
                                -out ${uid}_corrected.mzML
  """
}

/*
 * peptide identification
 * database search with MS-GF+ https://github.com/MSGFPlus/msgfplus
 */
process peptide_identification {

  tag "$uid"

  input:
    tuple uid, path(file_mzML) from ch_6
    path(file_fasta) from ch_database
    path(file_jar) from ch_msgfplus

  output:
   tuple uid, path("${uid}.idXML") into ch_7

  script:
  """
  MSGFPlusAdapter -in ${file_mzML} \\
                  -out ${uid}.idXML \\
                  -database ${file_fasta} \\
                  -executable ${file_jar} \\
                  -instrument Q_Exactive \\
                  -protocol none \\
                  -min_precursor_charge 2 \\
                  -max_precursor_charge 5 \\
                  -max_missed_cleavages 3 \\
                  -fixed_modifications 'Carbamidomethyl (C)' \\
                  -variable_modifications 'Oxidation (M)' 'Phospho (S)' \\
                  -java_memory 7000
  """
}

/*
 * peptide indexing
 * annotate PSMs with meta data such as target/decoy
 */
process peptide_indexing {

  tag "$uid"

  input:
    tuple uid, path(file_idXML) from ch_7
    path(file_fasta) from ch_database

  output:
    tuple uid, path("${uid}_indexed.idXML") into ch_8

  script:
  """
  PeptideIndexer -in ${file_idXML} \\
                 -out ${uid}_indexed.idXML \\
                 -fasta ${file_fasta} \\
                 -decoy_string DECOY_ \\
                 -missing_decoy_action warn \\
                 -write_protein_sequence \\
                 -write_protein_description \\
                 -enzyme:name Trypsin/P \\
                 -enzyme:specificity full
  """
}

ch_8.into { ch_9; ch_10 }

/*
 * false discovery rate
 * filter for false discovery rate
 */
process false_discovery_rate {

  tag "$uid"

  input:
    tuple uid, path(file_idXML) from ch_9

  output:
    tuple uid, path("${uid}_fdr.idXML") into ch_11

  script:
  """
  FalseDiscoveryRate -in ${file_idXML} \\
                     -out ${uid}_fdr.idXML \\
                     -PSM true \\
                     -FDR:PSM 0.01 \\
                     -protein false
  """
}

ch_11.into { ch_12; ch_13 }

/*
 * mzTab export of peptide identifications
 */
process mztab_export_peptide_id {

  tag "$uid"
  publishDir params.output, mode: 'copy'

  input:
    tuple uid, path(file_idXML) from ch_12

  output:
    tuple uid, path("${uid}_id.mzTab") into ch_14

  script:
  """
  MzTabExporter -in ${file_idXML} \\
                -out ${uid}_id.mzTab
  """
}

/*
 * map peptide IDs to peptide features
 */
process id_mapping {

  tag "$uid"

  input:
    tuple uid, path(file_featureXML), path(file_idXML) from ch_5.join(ch_13)

  output:
    path("${uid}_annotated.featureXML") into ch_15

  script:
  """
  IDMapper -in ${file_featureXML} \\
           -id ${file_idXML} \\
           -out ${uid}_annotated.featureXML \\
           -rt_tolerance 40.0 \\
           -mz_tolerance 10.0 \\
           -mz_measure ppm \\
           -mz_reference precursor \\
           -feature:use_centroid_rt true \\
           -feature:use_centroid_mz true
  """
}

/*
 * link peptide feature maps of all samples
 */
process link_maps {

  input:
    path files_featureXML from ch_15.collect()

  output:
    path "out.consensusXML" into ch_16

  script:
  """
  FeatureLinkerUnlabeledQT -in ${(files_featureXML as List).join(' ')}  \\
           -out out.consensusXML \\
           -algorithm:distance_RT:max_difference 200.0 \\
           -algorithm:distance_MZ:max_difference 20.0 \\
           -algorithm:distance_MZ:unit ppm
  """
}

/*
 * clean up IDs
 * The linked peptide features may contain different, contradicting IDs.
 * Keep only ID from highest intensity feature.
 */
process clean_up_ids {

  input:
    path file_consensusXML from ch_16

  output:
    path "out.consensusXML" into ch_17

  script:
  """
  IDConflictResolver -in ${file_consensusXML}  \\
           -out out.consensusXML \\
           -resolve_between_features highest_intensity
  """
}

/*
 * normalize peptide feature intensities in linked maps
 */
process normalize {

  input:
    path file_consensusXML from ch_17

  output:
    path "out.consensusXML" into ch_18

  script:
  """
  ConsensusMapNormalizer -in ${file_consensusXML}  \\
                         -out out.consensusXML \\
                         -algorithm_type quantile
  """
}

ch_18.into { ch_19; ch_20 }

/*
 * mzTab export of peptide quantifications
 */
process mztab_export_peptide_quant {

  publishDir params.output, mode: 'copy'

  input:
    path file_consensusXML from ch_20

  output:
    path "peptides_quant.mzTab" into ch_21

  script:
  """
  MzTabExporter -in ${file_consensusXML} \\
                -out peptides_quant.mzTab
  """
}

/*
 * ID posterior error probablility
 */
process error_probability {

  tag "$uid"

  input:
    tuple uid, path(file_idXML) from ch_10

  output:
    path("${uid}_err.idXML") into ch_22

  script:
  """
  IDPosteriorErrorProbability -in ${file_idXML} \\
                              -out ${uid}_err.idXML
  """
}

/*
 * merge IDs of all samples
 */
process merge_ids {

  input:
    path files_idXML from ch_22.collect()

  output:
    path "out.idXML" into ch_23

  script:
  """
  IDMerger -in ${(files_idXML as List).join(' ')}  \\
           -out out.idXML
  """
}

/*
 * protein inference
 */
process protein_inference {

  input:
    path file_idXML from ch_23

  output:
    path "out_inf.idXML" into ch_24

  script:
  """
  Epifany -in ${file_idXML}  \\
          -out out_inf.idXML
  """
}

/*
 * protein quantification
 */
process protein_quantification {

  publishDir params.output, mode: 'copy'

  input:
    path file_consensusXML from ch_19
    path file_idXML from ch_24

  output:
    path "out.csv" into ch_25

  script:
  """
  ProteinQuantifier -in ${file_idXML}  \\
                    -out out.csv \\
                    -mztab out.mzTab
  """
}
