# NOISY effectivity analysis 
This repository was created for a project testing the effectiveness of the NOISY program. 

You may find 3 categories of files here:
  * Python scripts and notebooks
  * Directories
  * Rest of files

## 1. Python scripts :snake:
  * __swisspfam_to_db_new.py__ — function for the separation of a swisspfam database raw text file into multiple tables in mariaDB instance
  * __pfam_table_creation.py__ — function for the creation of a table with mappings between tax_ids and mnemonics based on the speclist.txt file
  * __noisy_ds_construction_functions.py__ — functions that are necessary for the correct construction of the datasets (used by __/ds_constuctors/__ scripts)
  * __noisy_output_analysis_functions.py__ — functions that are required for the analysis of a snakemake pipeline output
  * __noisy_reference_tree_construction_functions.py__ — functions that are required for the correct reference tree construction
  * __noisy_tree_construction.ipynb__ — notebook that allow to construct reference trees based on a given Taxonomy dataframe
  * __noisy_statistics_and_calculations_new.ipynb__ — notebook that allow to calculate statistics for noisy runs and visualize it

## 2. Directories :file_folder:
  * __RF_statistics_df__ — folder with 6 .csv files with aggregated results for original/denoised runs and precalculated statistics
  * __ds_constructors__ — folder with 6 .py files aimed for the construction of datasets based on given MariaDB Taxonomy table and input parameters
  * __mapping_tables__ — folder with 6 .csv tables that contain aggregated information about domains, species mnemonic and tax_id of an aggregation class for every run from ds_constructors
  * __reference_trees__ — folder with 6 .nwk reference trees for species of each major subgroup, created with __noisy_tree_construction.ipynb__ notebook
  * __snakemake_configs__ — folder with 12 .yaml configs for every snakemake run (denoised and original for every group, bootstrap runs use the same config)

## 3. Rest of files :woozy_face:
  * __noisy_env.yml__ — .yml environment used for the project
  * __Snakefile__ — snakefile used for testing pipeline runs (configfile from __snakemake_configs__ and type of run (_denoised_ or _original_) may be specified)
  * __Snakefile_bootstrap__ — same snakefile but with additional option for bootstrap runs
  * __pipeline_runner.sh__ — file that allow to execute all of the snakemake pipeline runs at the same time
  * __pipeline_runner_bootstrap.sh__ — same file, but for bootstrap directories

## Additional comments and side notes :warning:
  * Please note, that on the current stage the project is not 100% reproducible on one click. Additional tuning of the relative and absolute paths will be required. 
  * TODO: Project scheme and description.
