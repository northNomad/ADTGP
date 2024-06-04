#' Raw protein data from MissionBio's 2PBMC dataset
#'
#' A dataset containing raw counts of 45 proteins from 3,751 cells
#' (a mix of two PBMC samples).
#' Sequencing was performed on an Illumina NextSeq 2000 instrument.
#' Protein reads were generated at 339x coverage of reads/antibody/cell.
#'
#' @format A matrix with 45 columns and 3,751 rows.
#' @source \url{https://portal.missionbio.com/datasets/2-PBMC-mix}
"missionbio_pbmc"

#' Raw protein data from isozyme switching dataset
#'
#' A dataset containing raw counts of 46 proteins from 979 cells.
#' An IDH1-mutated AML sample was mixed with an IDH2-mutated AML sample
#' at a 1-to-1 ratio. The mixed cells were injected into the tail vein
#'  of immunodeficient mice to generate mixed patient derived xenografts.
#'  The animals were treated with vehicle control or ivosidenib,
#'  an IDH1 inhibitor, to induce the differentiation of AML cells for eight
#'  weeks. At endpoint, bone marrow cells were harvested from the
#'  animals (n=5 per treatment) and processed through MissionBio's
#'  Tapestri single-cell multi-omics pipeline.
#'
#'  The "Patient" column was determined using the genotype of the cells, which
#'  is not provided in this dataset.
#'  veh and ivo are short-hand for vehicle and ivosidenib, respectively.
#'
#' @format A dataframe with 48 columns and 979 rows.
"idh_pdx"
