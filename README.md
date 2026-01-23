# GeneLab Metagenomics Sequencing Data Processing Workflow

> GeneLab, part of NASA's [Open Science Data Repository (OSDR)](https://www.nasa.gov/osdr/), has wrapped each step of the Low biomass [long-read](https://github.com/nasa/GeneLab_Data_Processing/blob/DEV_Metagenomics_low_biomass/Metagenomics/Low_Biomass/Nanopore) and [short-read](https://github.com/nasa/GeneLab_Data_Processing/blob/DEV_Metagenomics_low_biomass/Metagenomics/Low_Biomass/Illumina) metagenomics sequencing data processing pipelines into a Nextflow workflow. This repository contains information about the workflow along with instructions for installation and usage. Exact workflow run info and pipeline version used to process specific datasets hosted on the [OSDR data repository](https://osdr.nasa.gov/bio/repo/) are provided alongside their processed data in OSDR under 'Files' -> 'GeneLab Processed Metagenomics Files' -> 'Processing Info'. 

<br>

<p align="center">
<a href="images/GL-amplicon-subwayplot.pdf"><img src="images/GL-amplicon-subwayplot.png"></a>
</p>

<br>

## General Workflow Info

### Implementation Tools

The current GeneLab metagenomics sequencing data processing pipelines, Low biomass [long-read](https://github.com/nasa/GeneLab_Data_Processing/blob/DEV_Metagenomics_low_biomass/Metagenomics/Low_Biomass/Nanopore) and Low biomass [short-read](https://github.com/nasa/GeneLab_Data_Processing/blob/DEV_Metagenomics_low_biomass/Metagenomics/Low_Biomass/Illumina), are implemented as a [Nextflow](https://nextflow.io/) DSL2 workflow and utilizes [Singularity](https://docs.sylabs.io/guides/3.10/user-guide/introduction.html) containers, [Docker](https://docs.docker.com/get-started/) containers, or [conda](https://docs.conda.io/en/latest/) environments to install/run all tools. This workflow is run using the command line interface (CLI) of any unix-based system.  While knowledge of creating workflows in Nextflow is not required to run the workflow as is, [the Nextflow documentation](https://nextflow.io/docs/latest/index.html) is a useful resource for users who want to modify and/or extend this workflow.   

### Resource Requirements <!-- omit in toc -->

The table below details the default maximum resource allocations for individual Nextflow processes.

| CPU Cores | Memory |
|--------------------|------------------|
| 10                 | 600 GB           |

> **Note:** These per-process resource allocations are defaults. They can be adjusted by modifying `cpus` and `memory` directives in the configuration file: [`nextflow.config`](nextflow.config).

<br>

## Utilizing the Workflow

1. [Install Nextflow, Singularity, and Conda](#1-install-nextflow-singularity-and-conda)  
   1a. [Install Nextflow and Conda](#1a-install-nextflow-and-conda)  
   1b. [Install Singularity](#1b-install-singularity)  

2. [Download the Workflow Files](#2-download-the-workflow-files)  

3. [Fetch Singularity Images](#3-fetch-singularity-images)  

4. [Run the Workflows](#4-run-the-workflows)  
   4a. [Low Biomass Long Read Workflow](#4a-low-biomass-long-read-workflow)  
       4ai. [Approach 1: Approach 1: Start with pod5 or fast5 files as input](#4ai-approach-1-start-with-pod5-or-fast5-files-as-input)  
       4aii. [Approach 2: Start with multiple FASTQ files per sample as input](#4aii-approach-2-start-with-multiple-fastq-files-per-sample-as-input)  
       4aiii. [Approach 3: Start with one FASTQ file per sample as input](#4aiii-approach-3-start-with-one-fastq-file-per-sample-as-input)  
   4b. [Low Biomass Short Read Workflow](#4b-low-biomass-short-read-workflow)  
       4bi. [Approach 1: Start with paired-end FASTQ files as input](#4bi-approach-1-start-with-paired-end-fastq-files-as-input)  
       4bii. [Approach 2: Start with single-end FASTQ files as input](#4bii-approach-2-start-with-single-end-fastq-files-as-input)  
   4c. [Modify parameters and compute resources in the Nextflow config file](#4c-modify-parameters-and-compute-resources-in-the-nextflow-config-file)  

5. [Workflow Outputs](#5-workflow-outputs)  
   5a. [Main outputs](#5a-main-outputs)  
   5b. [Resource logs](#5b-resource-logs)  

<br>

---

### 1. Install Nextflow, Singularity, and Conda

#### 1a. Install Nextflow and Conda

Nextflow can be installed either through the [Anaconda bioconda channel](https://anaconda.org/bioconda/nextflow) or as documented on the [Nextflow documentation page](https://www.nextflow.io/docs/latest/getstarted.html).

> Note: If you wish to install conda, we recommend installing a Miniforge version appropriate for your system, as documented on the [conda-forge website](https://conda-forge.org/download/), where you can find basic binaries for most systems. More detailed miniforge documentation is available in the [miniforge github repository](https://github.com/conda-forge/miniforge).
> 
> Once conda is installed on your system, you can install the latest version of Nextflow by running the following commands:
> 
> ```bash
> conda install -c bioconda nextflow
> nextflow self-update
> ```
> You may also install [mamba](https://mamba.readthedocs.io/en/latest/index.html) first which is a faster implementation of conda and can be used as a drop-in replacement:
> ```bash
> conda install -c conda-forge mamba
> ```

<br>

#### 1b. Install Singularity

Singularity is a container platform that allows usage of containerized software. This enables the GeneLab workflow to retrieve and use all software required for processing without the need to install the software directly on the user's system.

We recommend installing Singularity on a system wide level as per the associated [documentation](https://docs.sylabs.io/guides/3.10/admin-guide/admin_quickstart.html).

> Note: Singularity is also available through the [Anaconda conda-forge channel](https://anaconda.org/conda-forge/singularity).

> Note: Alternatively, Docker can be used in place of Singularity. To get started with Docker, see the [Docker CE installation documentation](https://docs.docker.com/engine/install/).

<br>

---

### 2. Download the Workflow Files

All files required for utilizing the GeneLab metagenomics workflow for processing either short- or long-read low biomass sequencing data are available in this repository. To get a copy of the latest workflow version on to your system, clone this repository then `cd` into the repository directory by running the following commands: 

```bash
git clone https://github.com/olabiyi/NF-MetagenomeSeq.git
cd NF-MetagenomeSeq
```

<br>

---

### 3. Fetch Singularity Images

Although Nextflow can fetch Singularity images from a url, doing so may cause issues as detailed [here](https://github.com/nextflow-io/nextflow/issues/1210).

To avoid this issue, run the following command to fetch the Singularity images prior to running the GeneLab metagenomics workflow:

> Note: This command should be run from within the `NF-MetagenomeSeq` directory that was downloaded in [step 2](#2-download-the-workflow-files) above. Depending on your network speed, fetching the images will take ~20 minutes. Approximately 4GB of RAM is needed to download and build the Singularity images.

```bash
bash ./bin/prepull_singularity.sh nextflow.config
```

Once complete, a `singularity` folder containing the Singularity images will be created. Run the following command to export this folder as a Nextflow configuration environment variable to ensure Nextflow can locate the fetched images:

```bash
export NXF_SINGULARITY_CACHEDIR=$(pwd)/singularity
```

<br>

---

### 4. Run the Workflows

> **Notes:**
> - All the commands in this step assume that the workflow will be run from within the `NF-MetagenomeSeq` directory that was downloaded in [step 2](#2-download-the-workflow-files) above. They may also be run from a different location by providing the full paths to the launch.sh script, low_biomass_nanopore.nf, low_biomass_illumina.nf, and nextflow.config workflow files in the `NF-MetagenomeSeq` directory.*
>
> - Nextflow commands use both single hyphen arguments (e.g. -help) that denote general Nextflow arguments and double hyphen arguments (e.g. --input_file) that denote workflow specific parameters.  Take care to use the proper number of hyphens for each argument.

<br>

#### 4a. Low Biomass Long Read Workflow

The GeneLab Metagenomics Low Biomass Long Read workflow is designed to process data generated from long-read platforms such as [Oxford Nanopore](https://nanoporetech.com/) using the [GeneLab Metagenomics Low Biomass Long Read Pipeline](https://github.com/nasa/GeneLab_Data_Processing/blob/DEV_Metagenomics_low_biomass/Metagenomics/Low_Biomass/Nanopore). Below are 3 different approaches for running the workflow, depending on the input files provided.

For options and detailed help on how to run the low biomass long read workflow, run the following command:

```bash
nextflow run low_biomass_nanopore.nf --help
```

##### 4ai. Approach 1: Start with pod5 or fast5 files as input

```bash
bash ./launch.sh processing \
   low_biomass_nanopore.nf \
   '--input_file input_dir_barcodes.csv --errorStrategy "ignore" '
```

<br>

##### 4aii. Approach 2: Start with multiple FASTQ files per sample as input

```bash
bash ./launch.sh processing \
   low_biomass_nanopore.nf \
   '--input_file multiple.csv --errorStrategy "ignore" '
```

<br>

##### 4aiii. Approach 3: Start with one FASTQ file per sample as input

```bash
bash ./launch.sh processing \
   low_biomass_nanopore.nf \
   '--input_file single.csv --errorStrategy "ignore" '
```

<br>

**Required Parameters For All Long Read Approaches:**

* `./launch.sh` - Positional argument specifying the bash file that contains the nextflow command to execute. If running in a directory other than `NF-MetagenomeSeq`, replace with the full path to the launch.sh script.
* `processing` - Positional argument specifying the run mode, `processing` refers to the mode used to process the input data using the respective GeneLab pipeline.
* `low_biomass_nanopore.nf` - Instructs Nextflow to run the Genelab Low Biomass Long Read (Nanopore) workflow. If running in a directory other than `NF-MetagenomeSeq`, replace with the full path to the low_biomass_nanopore.nf workflow file.
* `'--input_file *.csv --errorStrategy "ignore" '` - Positional argument specifying a set of parameters to pass into the nextflow run command.
  * `--errorStrategy "ignore"` - Instructions nextflow to continue processing the dataset even if an error is encountered.
  * `--input_file *.csv` - Specifies the input csv file containing required metadata about the samples including barcode information and paths to the input file(s) for each sample.
  * > *Note: This input files requires specifc formatting to be interpreted correctly. Please see the [runsheet documentation](examples/runsheet) in this repository for examples on how to format this file type for each approach.

<br>

#### 4b. Low Biomass Short Read Workflow

The GeneLab Metagenomics Low Biomass Short Read workflow is designed to process data generated from short-read platforms such as [Illumina](https://www.illumina.com/) using the [GeneLab Metagenomics Low Biomass Short Read Pipeline](https://github.com/nasa/GeneLab_Data_Processing/blob/DEV_Metagenomics_low_biomass/Metagenomics/Low_Biomass/Illumina). Below are 2 different approaches for running the workflow, depending on the input files provided.

For options and detailed help on how to run the low biomass short read workflow, run the following command:

```bash
nextflow run low_biomass_illumina.nf --help
```

##### 4bi. Approach 1: Start with paired-end FASTQ files as input

```bash
bash ./launch.sh processing \
   low_biomass_illumina.nf \
   '--input_file PE_file.csv --errorStrategy "ignore" --technology "illumina" '
```

<br>

##### 4bii. Approach 2: Start with single-end FASTQ files as input

```bash
bash ./launch.sh processing \
   low_biomass_illumina.nf \
   '--input_file SE_file.csv --errorStrategy "ignore" --technology "illumina" '
```

<br>

**Required Parameters For All Short Read Approaches:**

* `./launch.sh` - Positional argument specifying the bash file that contains the nextflow command to execute. If running in a directory other than `NF-MetagenomeSeq`, replace with the full path to the launch.sh script.
* `processing` - Positional argument specifying the run mode, `processing` refers to the mode used to process the input data using the respective GeneLab pipeline.
* `low_biomass_illumina.nf` - Instructs Nextflow to run the Genelab Low Biomass Short Read (Illumina) workflow. If running in a directory other than `NF-MetagenomeSeq`, replace with the full path to the low_biomass_illumina.nf workflow file.
* `'--input_file *.csv --errorStrategy "ignore" --technology "illumina" '` - Positional argument specifying a set of parameters to pass into the nextflow run command.
  * `--errorStrategy "ignore"` - Instructions nextflow to continue processing the dataset even if an error is encountered.
  * `--technology "illumina"` - Specifies the technology type used to generate the sequencing data.
  * `--input_file *.csv` - Specifies the input csv file containing required metadata about the samples including paths to the input file(s) for each sample.
  * > *Note: This input files requires specifc formatting to be interpreted correctly. Please see the [runsheet documentation](examples/runsheet) in this repository for examples on how to format this file type for each approach.

<br>

**Additional [Optional] Parameters For All Approaches For Both Long- and Short-Read**
> *Note: See `nextflow run -h` and [Nextflow's CLI run command documentation](https://nextflow.io/docs/latest/cli.html#run) for more options and details on how to run Nextflow.*

* `'-profile'` – Specifies the configuration profile(s) to load (multiple options can be provided as a comma-separated list).
   * Software environment profile options (choose one):
      * `singularity` - instructs Nextflow to use Singularity container environments
      * `docker` - instructs Nextflow to use Docker container environments
      * `conda` - instructs Nextflow to use conda environments via the conda package manager
        > *Note: By default, Nextflow will create environments at runtime using the `singularity` approach.*
* `'--assay_suffix ""'` – Specifies the suffix to add to each output file.

<br>

#### 4c. Modify parameters and compute resources in the Nextflow config file

Additionally, all parameters and workflow resources can be directly specified in the [nextflow.config](./workflow_code/nextflow.config) file. For detailed instructions on how to modify and set parameters in the config file, please see the [documentation here](https://www.nextflow.io/docs/latest/config.html).

Once you've downloaded the workflow template, you can modify the parameters in the `params` scope and cpus/memory requirements in the `process` scope in your downloaded version of the [nextflow.config](workflow_code/nextflow.config) file as needed in order to match your dataset and system setup. Additionally, if necessary, you can modify each variable in the [nextflow.config](workflow_code/nextflow.config) file to be consistent with the study you want to process and the computer you're using for processing.

<br>

---

### 5. Workflow Outputs

#### 5a. Main Outputs

* The outputs from the GeneLab Low Biomass Long Read Metagenomics workflow are documented in the [GL-DPPD-7116](https://github.com/nasa/GeneLab_Data_Processing/blob/DEV_Metagenomics_low_biomass/Metagenomics/Low_Biomass/Nanopore/GL-DPPD-7116.md) processing pipeline.
  
* The outputs from the GeneLab Low Biomass Short Read Metagenomics workflow are documented in the [GL-DPPD-7117](https://github.com/nasa/GeneLab_Data_Processing/blob/DEV_Metagenomics_low_biomass/Metagenomics/Low_Biomass/Illumina/GL-DPPD-7117.md) processing pipeline.

#### 5b. Resource Logs

Standard Nextflow resource usage logs are also produced as follows:

**Nextflow Resource Usage Logs**
   - Resource_Usage/execution_report_{timestamp}.html (an html report that includes metrics about the workflow execution including computational resources and exact workflow process commands)
   - Resource_Usage/execution_timeline_{timestamp}.html (an html timeline for all processes executed in the workflow)
   - Resource_Usage/execution_trace_{timestamp}.txt (an execution tracing file that contains information about each process executed in the workflow, including: submission time, start time, completion time, cpu and memory used, machine-readable output)

> Further details about these logs can also found within [this Nextflow documentation page](https://www.nextflow.io/docs/latest/tracing.html#execution-report).

<br>

---

## License

The software for the GeneLab Metagenomics workflow is released under the [NASA Open Source Agreement (NOSA) Version 1.3](License/Metagenomics_NOSA_License.pdf).


### 3rd Party Software Licenses

Licenses for the 3rd party open source software utilized in the GeneLab Metagenomics workflow can be found in the [3rd_Party_Licenses sub-directory](License/3rd_Party_Licenses). 

<br>

---

## Notices

Copyright © 2021 United States Government as represented by the Administrator of the National Aeronautics and Space Administration.  All Rights Reserved.

### Disclaimers

No Warranty: THE SUBJECT SOFTWARE IS PROVIDED "AS IS" WITHOUT ANY WARRANTY OF ANY KIND, EITHER EXPRESSED, IMPLIED, OR STATUTORY, INCLUDING, BUT NOT LIMITED TO, ANY WARRANTY THAT THE SUBJECT SOFTWARE WILL CONFORM TO SPECIFICATIONS, ANY IMPLIED WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE, OR FREEDOM FROM INFRINGEMENT, ANY WARRANTY THAT THE SUBJECT SOFTWARE WILL BE ERROR FREE, OR ANY WARRANTY THAT DOCUMENTATION, IF PROVIDED, WILL CONFORM TO THE SUBJECT SOFTWARE. THIS AGREEMENT DOES NOT, IN ANY MANNER, CONSTITUTE AN ENDORSEMENT BY GOVERNMENT AGENCY OR ANY PRIOR RECIPIENT OF ANY RESULTS, RESULTING DESIGNS, HARDWARE, SOFTWARE PRODUCTS OR ANY OTHER APPLICATIONS RESULTING FROM USE OF THE SUBJECT SOFTWARE.  FURTHER, GOVERNMENT AGENCY DISCLAIMS ALL WARRANTIES AND LIABILITIES REGARDING THIRD-PARTY SOFTWARE, IF PRESENT IN THE ORIGINAL SOFTWARE, AND DISTRIBUTES IT "AS IS."

Waiver and Indemnity:  RECIPIENT AGREES TO WAIVE ANY AND ALL CLAIMS AGAINST THE UNITED STATES GOVERNMENT, ITS CONTRACTORS AND SUBCONTRACTORS, AS WELL AS ANY PRIOR RECIPIENT.  IF RECIPIENT'S USE OF THE SUBJECT SOFTWARE RESULTS IN ANY LIABILITIES, DEMANDS, DAMAGES, EXPENSES OR LOSSES ARISING FROM SUCH USE, INCLUDING ANY DAMAGES FROM PRODUCTS BASED ON, OR RESULTING FROM, RECIPIENT'S USE OF THE SUBJECT SOFTWARE, RECIPIENT SHALL INDEMNIFY AND HOLD HARMLESS THE UNITED STATES GOVERNMENT, ITS CONTRACTORS AND SUBCONTRACTORS, AS WELL AS ANY PRIOR RECIPIENT, TO THE EXTENT PERMITTED BY LAW.  RECIPIENT'S SOLE REMEDY FOR ANY SUCH MATTER SHALL BE THE IMMEDIATE, UNILATERAL TERMINATION OF THIS AGREEMENT.
