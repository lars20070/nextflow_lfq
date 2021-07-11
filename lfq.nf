#!/usr/bin/env nextflow

// directory with mzML files
params.input = "$baseDir/data/*.mzML"

params.database = "$baseDir/data/W82_soybase_a2v1_and_pMOZ52.fasta"
params.msgfplus = "/Volumes/GoogleDrive/Shared drives/LCMS/Analyses/soybean analysis 20210527/MSGFPlus.jar"

/*
 * channels for peptide identification and peptide quantification
 * The channel contains pairs of an unique ID (base name of the file) and the mzML file itself.
 * Note that the unique ID is passed through all processes of the workflow in order to match corresponding files.
 */
Channel.fromPath(params.input).map { file -> [file.baseName, file ]}.into { ch_1; ch_2 }

/*
 * detect peptide features in MS1 spaectral data
 */
process peptide_feature_detection {

  //publishDir "$baseDir/data/results", mode: 'copy'

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
 * In that case, we shift the precursor position back to the mono-sicotopic peak.
 * This correction simplifies the subsequent database search.
 */
process precursor_correction {

  //publishDir "$baseDir/data/results", mode: 'copy'

  input:
    tuple uid, path(file_mzML), path(file_featureXML) from ch_2.join(ch_4)

  output:
    tuple uid, path("${uid}.mzML") into ch_6

  script:
  """
  HighResPrecursorMassCorrector -in ${file_mzML} \\
                                -feature:in ${file_featureXML} \\
                                -out ${uid}.mzML
  """
}

/*
 * peptide identification
 */
process peptide_identification {

  
}