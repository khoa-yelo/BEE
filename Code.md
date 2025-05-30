# 🧑‍💻 Code in Bioinformatics Engineering

Welcome to the coding section of **Bioinformatics Engineering for Everyone**! This module will give you hands-on experience with essential tools for modern bioinformatics.

## 📖 Overview

Modern bioinformatics requires managing complex codebases and running large-scale analyses. We'll learn through hands-on tutorials:
1. Version control with Git
2. High-performance computing with SLURM
3. Workflow management with Nextflow

## 🎯 Learning Goals
- Master Git for code version control and collaboration
- Run parallel jobs efficiently on clusters using SLURM arrays
- Build reproducible pipelines with Nextflow
- Apply best practices for bioinformatics code organization

## 🛠️ Tutorials

### 1. Git Essentials

**Why Use Git?**
- **Version Control**: Track changes in your code over time
  - Recover previous versions if something breaks
  - Understand when and why changes were made
  - Experiment safely with new features
- **Collaboration**: Work effectively with others
  - Multiple people can work on the same code
  - Track who made what changes
  - Resolve conflicts when changes overlap
- **Reproducibility**: Essential for scientific work
  - Record exact code version used for results
  - Share code with other researchers
  - Document development history
- **Common Scenarios**:
  - Developing analysis pipelines
  - Writing research scripts
  - Collaborating on code with lab members
  - Publishing code with papers

#### Basic Git Workflow
```bash
# Initialize a repository
git init

# Stage changes
git add <filename>
git add .  # add all changes

# Commit changes
git commit -m "descriptive message"

# Check status
git status
```

#### Branching and Collaboration
```bash
# Create and switch to a new branch
git checkout -b feature_name

# Push to remote
git push origin feature_name

# Pull updates
git pull origin main
```

#### Best Practices
- Write clear commit messages
- Use branches for new features
- Create .gitignore for data files
- Regular commits, regular pushes

### 2. SLURM Array Jobs

**Why Use SLURM Arrays?**
- **Parallel Processing**: Handle large-scale analyses
  - Process hundreds/thousands of samples
  - Utilize cluster resources efficiently
  - Reduce total processing time
- **Resource Management**: Optimize compute usage
  - Distribute work across nodes
  - Control memory and CPU allocation
  - Schedule jobs efficiently
- **Workflow Control**: Manage complex analyses
  - Track job progress
  - Handle failures gracefully
  - Organize output files
- **Common Scenarios**:
  - Processing multiple sequencing samples
  - Running analyses with different parameters
  - Batch processing of large datasets
  - Parallel simulations or analyses

#### Example: Processing FastQ Files
```bash
#!/bin/bash

#SBATCH --job-name=fastq_proc
#SBATCH --array=0-3          # 4 parallel jobs (0,1,2,3)
#SBATCH --output=logs/%A_%a.log
#SBATCH --error=logs/%A_%a.err
#SBATCH --mem=4G
#SBATCH --cpus-per-task=1

# Exit on error
set -e

# Setup directories
FASTQ_DIR="data/fastq"
OUT_DIR="results/processed"
mkdir -p ${OUT_DIR} logs

# List all fastq files and get total count
files=(${FASTQ_DIR}/*.fastq.gz)
num_files=${#files[@]}
num_array_jobs=4  # Same as array range (0-3)

echo "Total files: ${num_files}"
echo "This is job ${SLURM_ARRAY_TASK_ID} out of ${num_array_jobs}"

# Process files distributed by modulo operation
for ((i=0; i<${num_files}; i++)); do
    # Only process files where index modulo matches this job's array ID
    if [ $((i % num_array_jobs)) -eq ${SLURM_ARRAY_TASK_ID} ]; then
        fastq=${files[$i]}
        basename=$(basename ${fastq} .fastq.gz)
        
        echo "[Job ${SLURM_ARRAY_TASK_ID}] Processing ${basename}"
        
        # Example: Run FastQC on each file
        fastqc \
            --outdir ${OUT_DIR} \
            --threads ${SLURM_CPUS_PER_TASK} \
            ${fastq}
            
        echo "[Job ${SLURM_ARRAY_TASK_ID}] Completed ${basename}"
    fi
done

echo "All files for job ${SLURM_ARRAY_TASK_ID} completed"
```

#### How It Works

1. **Job Configuration**
   - Creates 4 parallel jobs (array indices 0-3)
   - Each job gets 4GB memory and 1 CPU
   - Separate log files for each array task

2. **File Distribution**
   - Lists all `.fastq.gz` files in input directory
   - Distributes files across jobs using modulo
   - If you have 10 files:
     - Job 0 processes files 0,4,8
     - Job 1 processes files 1,5,9
     - Job 2 processes files 2,6
     - Job 3 processes files 3,7

3. **Error Handling**
   - `set -e` exits on any error
   - Creates output directories
   - Logs progress for each file

#### Common Commands
```bash
# Submit job
sbatch process_fastq.sh

# Check job status
squeue -u $USER

# Cancel job
scancel <job_id>

# Monitor progress
tail -f logs/<job_id>_<array_id>.log
```

#### Best Practices
- Test script with small dataset first
- Monitor resource usage
- Use appropriate array size for your data
- Include error handling and logging
- Create output directories before processing

### 3. Nextflow Pipelines

**Why Use Nextflow?**
- **Reproducibility**: Ensure consistent analyses
  - Capture entire workflow
  - Include software dependencies
  - Document all parameters
- **Scalability**: Handle complex workflows
  - Run locally or on clusters
  - Process thousands of samples
  - Optimize resource usage
- **Portability**: Work anywhere
  - Run same pipeline on different systems
  - Container integration
  - Cloud compatibility
- **Common Scenarios**:
  - RNA-seq analysis pipelines
  - Variant calling workflows
  - Multi-step data processing
  - Integration of multiple tools

#### Setting Up a Nextflow Project

1. **Project Structure**
```bash
my_pipeline/
├── main.nf           # Main pipeline script
├── nextflow.config   # Pipeline configuration
├── params.yml        # Parameter definitions
├── modules/          # Reusable process definitions
│   ├── fastqc.nf
│   ├── trim.nf
│   └── align.nf
├── conf/            # Configuration profiles
│   ├── base.config
│   ├── slurm.config
│   └── docker.config
├── bin/             # Binary scripts/tools
├── lib/             # Custom functions
└── assets/          # Additional resources
```

2. **Configuration Files**

**nextflow.config**
```groovy
// Load parameter file
params.config = 'params.yml'

// Include configuration profiles
profiles {
    standard {
        process.executor = 'local'
    }
    slurm {
        process.executor = 'slurm'
        includeConfig 'conf/slurm.config'
    }
    docker {
        docker.enabled = true
        includeConfig 'conf/docker.config'
    }
}

// Default process configuration
process {
    cpus = 1
    memory = '2 GB'
    errorStrategy = 'retry'
    maxRetries = 3
}
```

**params.yml**
```yaml
# Input/Output
input_dir: "data/fastq"
output_dir: "results"

# Reference files
genome: "ref/genome.fa"
annotation: "ref/genes.gtf"

# Tool parameters
min_quality: 20
min_length: 50
threads: 4
```

**conf/slurm.config**
```groovy
process {
    executor = 'slurm'
    queue = 'normal'
    
    withLabel: small {
        cpus = 1
        memory = '4 GB'
        time = '1h'
    }
    
    withLabel: medium {
        cpus = 4
        memory = '16 GB'
        time = '6h'
    }
    
    withLabel: large {
        cpus = 8
        memory = '32 GB'
        time = '12h'
    }
}
```

#### Key Nextflow Concepts

1. **Channels**
```groovy
// Create channels from files
reads_ch = Channel
    .fromFilePairs("data/*_{1,2}.fastq.gz")
    .map { id, files -> tuple(id, files[0], files[1]) }

// Value channels
params_ch = Channel.value(params.min_quality)

// Combine channels
combined_ch = reads_ch.combine(params_ch)
```

2. **Processes**
```groovy
// Basic process structure
process FASTQC {
    label 'small'
    tag "$sample_id"
    
    input:
    tuple val(sample_id), path(read1), path(read2)
    
    output:
    tuple val(sample_id), path("*_fastqc.html")
    
    script:
    """
    fastqc $read1 $read2
    """
}
```

3. **Workflows**
```groovy
// Workflow composition
workflow {
    FASTQC(reads_ch)
    TRIM(reads_ch)
    ALIGN(TRIM.out)
}

// Conditional workflows
workflow {
    if (params.run_qc) {
        FASTQC(reads_ch)
    }
    ALIGN(reads_ch)
}
```

#### Complete Example: RNA-seq Pipeline

```groovy
#!/usr/bin/env nextflow

// Parameter defaults
params.reads = "data/*_{1,2}.fastq.gz"
params.genome = "ref/genome.fa"
params.annotation = "ref/genes.gtf"
params.outdir = "results"

// Log parameters
log.info """\
    RNA-SEQ PIPELINE
    ================
    reads: ${params.reads}
    genome: ${params.genome}
    annotation: ${params.annotation}
    outdir: ${params.outdir}
    """

// Input channels
reads_ch = Channel
    .fromFilePairs(params.reads, checkIfExists: true)
    .map { id, files -> tuple(id, files[0], files[1]) }

// Quality Control
process FASTQC {
    label 'small'
    publishDir "${params.outdir}/fastqc", mode: 'copy'
    
    input:
    tuple val(sample_id), path(read1), path(read2)
    
    output:
    tuple val(sample_id), path("*_fastqc.html")
    
    script:
    """
    fastqc -t ${task.cpus} $read1 $read2
    """
}

// Trimming
process TRIM {
    label 'medium'
    publishDir "${params.outdir}/trimmed", mode: 'copy'
    
    input:
    tuple val(sample_id), path(read1), path(read2)
    
    output:
    tuple val(sample_id), path("*_trimmed_1.fq.gz"), path("*_trimmed_2.fq.gz")
    
    script:
    """
    trimmomatic PE -threads ${task.cpus} \
        $read1 $read2 \
        ${sample_id}_trimmed_1.fq.gz ${sample_id}_unpaired_1.fq.gz \
        ${sample_id}_trimmed_2.fq.gz ${sample_id}_unpaired_2.fq.gz \
        LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36
    """
}

// Alignment
process ALIGN {
    label 'large'
    publishDir "${params.outdir}/aligned", mode: 'copy'
    
    input:
    tuple val(sample_id), path(trimmed1), path(trimmed2)
    path genome
    
    output:
    tuple val(sample_id), path("*.bam"), path("*.bai")
    
    script:
    """
    STAR --runThreadN ${task.cpus} \
        --genomeDir $genome \
        --readFilesIn $trimmed1 $trimmed2 \
        --outFileNamePrefix ${sample_id}_ \
        --outSAMtype BAM SortedByCoordinate
    
    samtools index ${sample_id}_Aligned.sortedByCoord.out.bam
    """
}

// Main workflow
workflow {
    // Quality control
    FASTQC(reads_ch)
    
    // Trimming and alignment
    trimmed_reads = TRIM(reads_ch)
    aligned_reads = ALIGN(trimmed_reads, params.genome)
}
```

#### Running the Pipeline

```bash
# Run locally
nextflow run main.nf -profile standard

# Run on SLURM
nextflow run main.nf -profile slurm

# Resume failed run
nextflow run main.nf -resume

# Override parameters
nextflow run main.nf --reads "new_data/*_{1,2}.fastq.gz"
```

## 🔧 Practical Examples

### Example 1: RNA-seq Pipeline
```groovy
// RNA-seq workflow example
workflow {
    read_pairs_ch = Channel.fromFilePairs(params.reads)
    
    FASTQC(read_pairs_ch)
    TRIM(read_pairs_ch)
    ALIGN(TRIM.out)
    COUNT(ALIGN.out)
}
```

### Example 2: Variant Calling Array
```bash
#!/bin/bash
#SBATCH --array=1-22 # One job per chromosome

# Process one chromosome
chr=${SLURM_ARRAY_TASK_ID}
bcftools call -v -m -O z \
    -o results/chr${chr}.vcf.gz \
    input/chr${chr}.bcf
```

## 🎓 Practice Exercises

1. **Git Exercise**
   - Fork a repository
   - Make changes in a branch
   - Submit a pull request

2. **SLURM Exercise**
   - Write an array job script
   - Process multiple input files
   - Monitor resource usage

3. **Nextflow Exercise**
   - Create a basic pipeline
   - Add error handling
   - Use different profiles

Ready to code? Start with the Git tutorial and work your way through each section! 