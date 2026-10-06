# Input data

This repository contains no raw or processed expression data. Public inputs are fetched from the accession/source URLs in config/data_sources.json and the stage-specific configuration files. Preserve repository licences and source checksums.

The H9c2 input is private and is not supplied. No synthetic expression matrix is substituted. To run the private branch, an authorized user must place the original Summary.tar.gz and Report.tar.gz in a directory specified by H9C2_SOURCE_DIR or config/runtime.local.json. Module 08 lists exact required archive members; config/local_rat.json specifies six C and six L biological samples. H9c2-dependent modules cannot reproduce those results from this repository alone.

MSigDB Hallmark 2025.1.Hs must be obtained under the source terms and placed at the gmt path in config/biology_extension_v1.json. The gene sets are not redistributed here. The adult heart reference uses the public sources recorded in modules 93–97 and cellxgene_read_utils.py.
