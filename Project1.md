# 🦠 Project 1: Bacterial Genome Assembly and Analysis

## 🎯 Project Overview

In this project, you will:
1. Assemble a bacterial genome from long-read sequencing data
2. Annotate the genome to identify coding sequences
3. Cluster protein sequences to identify homologs
4. Analyze protein families with multiple copies

## 📋 Learning Objectives
- Master long-read genome assembly using Flye
- Understand bacterial genome annotation with Bakta
- Learn protein sequence clustering with MMseqs2
- Analyze and visualize protein families
- Work with Nextflow for pipeline development

## 🛠️ Tools Required
- Flye (v2.9.2+) for long-read assembly
- Bakta (v1.8+) for genome annotation
- MMseqs2 (v14.7e284+) for protein clustering
- Seqkit (v2.5+) for sequence manipulation
- R/Python for analysis and visualization

## 📊 Input Data
- PacBio/Oxford Nanopore long reads (FASTQ format)
- Bakta database
- Computing resources (16+ GB RAM recommended)

## 🔬 Project Pipeline

### Step 1: Setup Project Structure

```bash
project1/
├── data/
│   ├── raw_reads/
│   ├── bakta_db/
│   └── results/
├── scripts/
│   ├── plot_clusters.R
│   └── analyze_homologs.py
├── main.nf
├── nextflow.config
└── README.md
```

### Step 2: Nextflow Pipeline Implementation

Create `main.nf`:

```groovy
#!/usr/bin/env nextflow

// Parameters
params.reads = "data/raw_reads/*.fastq.gz"
params.outdir = "results"
params.bakta_db = "data/bakta_db"
params.min_contig_length = 1000
params.min_cluster_size = 2
params.min_seq_id = 0.9

// Log parameters
log.info """
    BACTERIAL GENOME ANALYSIS PIPELINE
    =================================
    reads        : ${params.reads}
    bakta_db    : ${params.bakta_db}
    outdir      : ${params.outdir}
    """

// Input channels
reads_ch = Channel
    .fromPath(params.reads)
    .map { file -> tuple(file.baseName, file) }

// 1. Genome Assembly
process ASSEMBLE {
    label 'large_mem'
    publishDir "${params.outdir}/assembly", mode: 'copy'
    
    input:
    tuple val(sample_id), path(reads)
    
    output:
    tuple val(sample_id), path("${sample_id}_assembly.fasta")
    
    script:
    """
    flye --nano-raw $reads \
        --out-dir assembly \
        --threads ${task.cpus} \
        --min-overlap 1000
    
    mv assembly/assembly.fasta ${sample_id}_assembly.fasta
    """
}

// 2. Genome Annotation
process ANNOTATE {
    label 'large_mem'
    publishDir "${params.outdir}/annotation", mode: 'copy'
    
    input:
    tuple val(sample_id), path(assembly)
    path bakta_db
    
    output:
    tuple val(sample_id), path("${sample_id}_proteins.faa"), path("${sample_id}_annotation.gff3")
    
    script:
    """
    bakta --db ${bakta_db} \
        --threads ${task.cpus} \
        --prefix ${sample_id} \
        --output annotation \
        ${assembly}
    
    mv annotation/${sample_id}.faa ${sample_id}_proteins.faa
    mv annotation/${sample_id}.gff3 ${sample_id}_annotation.gff3
    """
}

// 3. Protein Clustering
process CLUSTER {
    label 'medium_mem'
    publishDir "${params.outdir}/clusters", mode: 'copy'
    
    input:
    tuple val(sample_id), path(proteins), path(annotation)
    
    output:
    tuple val(sample_id), path("${sample_id}_clusters.tsv")
    
    script:
    """
    # Create MMseqs2 database
    mmseqs createdb ${proteins} proteins_db
    
    # Cluster sequences
    mmseqs cluster proteins_db cluster_db tmp \
        --min-seq-id ${params.min_seq_id} \
        --threads ${task.cpus}
    
    # Create TSV output
    mmseqs createtsv proteins_db proteins_db cluster_db ${sample_id}_clusters.tsv
    """
}

// 4. Analyze Homologs
process ANALYZE {
    label 'small_mem'
    publishDir "${params.outdir}/analysis", mode: 'copy'
    
    input:
    tuple val(sample_id), path(clusters)
    
    output:
    path "${sample_id}_homolog_analysis.pdf"
    path "${sample_id}_homolog_summary.tsv"
    
    script:
    """
    # Run analysis script
    analyze_homologs.py \
        --input ${clusters} \
        --min-size ${params.min_cluster_size} \
        --output ${sample_id}_homolog_summary.tsv
    
    # Generate visualization
    plot_clusters.R \
        --input ${sample_id}_homolog_summary.tsv \
        --output ${sample_id}_homolog_analysis.pdf
    """
}

// Main workflow
workflow {
    // Assembly
    assembled = ASSEMBLE(reads_ch)
    
    // Annotation
    annotated = ANNOTATE(assembled, params.bakta_db)
    
    // Clustering
    clustered = CLUSTER(annotated)
    
    // Analysis
    ANALYZE(clustered)
}
```

### Step 3: Configuration

Create `nextflow.config`:

```groovy
// Process resource configuration
process {
    // Default resources
    cpus = 1
    memory = '4 GB'
    
    withLabel: small_mem {
        cpus = 2
        memory = '8 GB'
        time = '2h'
    }
    
    withLabel: medium_mem {
        cpus = 4
        memory = '16 GB'
        time = '6h'
    }
    
    withLabel: large_mem {
        cpus = 8
        memory = '32 GB'
        time = '12h'
    }
}

// Profiles
profiles {
    standard {
        process.executor = 'local'
    }
    
    slurm {
        process.executor = 'slurm'
        process.queue = 'normal'
    }
}

// Docker containers
docker {
    enabled = true
    runOptions = '-u $(id -u):$(id -g)'
}

// Container configurations
process {
    withName: ASSEMBLE {
        container = 'nanozoo/flye:2.9.2'
    }
    withName: ANNOTATE {
        container = 'oschwengers/bakta:1.8.1'
    }
    withName: CLUSTER {
        container = 'soedinglab/mmseqs2:14.7e284'
    }
}
```

### Step 4: Analysis Scripts

Create `scripts/analyze_homologs.py`:

```python
#!/usr/bin/env python3

import pandas as pd
import argparse
from collections import Counter

def analyze_clusters(input_file, min_size, output_file):
    # Read clusters
    clusters = pd.read_csv(input_file, sep='\t', header=None)
    
    # Count proteins per cluster
    cluster_sizes = Counter(clusters[0])
    
    # Filter by minimum size
    large_clusters = {k:v for k,v in cluster_sizes.items() if v >= min_size}
    
    # Create summary
    summary = pd.DataFrame({
        'cluster_id': list(large_clusters.keys()),
        'size': list(large_clusters.values())
    })
    
    # Save results
    summary.to_csv(output_file, sep='\t', index=False)

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('--input', required=True)
    parser.add_argument('--min-size', type=int, default=2)
    parser.add_argument('--output', required=True)
    args = parser.parse_args()
    
    analyze_clusters(args.input, args.min_size, args.output)
```

Create `scripts/plot_clusters.R`:

```R
#!/usr/bin/env Rscript

library(ggplot2)
library(optparse)

# Parse arguments
option_list = list(
    make_option("--input", type="character"),
    make_option("--output", type="character")
)
opt = parse_args(OptionParser(option_list=option_list))

# Read data
data <- read.delim(opt$input)

# Create plot
p <- ggplot(data, aes(x=reorder(cluster_id, -size), y=size)) +
    geom_bar(stat="identity", fill="steelblue") +
    theme_minimal() +
    labs(
        title="Protein Family Sizes",
        x="Cluster ID",
        y="Number of Proteins"
    ) +
    theme(axis.text.x = element_text(angle=45, hjust=1))

# Save plot
ggsave(opt$output, p, width=10, height=6)
```

## 📝 Project Tasks

1. **Setup (10%)**
   - Clone repository
   - Install required tools
   - Download test dataset

2. **Pipeline Development (40%)**
   - Implement Nextflow pipeline
   - Configure resource requirements
   - Add error handling
   - Test with small dataset

3. **Analysis (30%)**
   - Run pipeline on full dataset
   - Analyze cluster sizes
   - Identify largest protein families
   - Visualize results

4. **Documentation (20%)**
   - Document pipeline usage
   - Explain parameter choices
   - Summarize findings
   - Suggest improvements

## 📊 Expected Outputs

1. Assembled genome (FASTA)
2. Genome annotation (GFF3)
3. Protein sequences (FAA)
4. Cluster assignments (TSV)
5. Analysis report (PDF)
6. Summary statistics (TSV)

## 💡 Tips for Success

1. **Quality Control**
   - Check assembly completeness
   - Validate protein predictions
   - Verify cluster assignments

2. **Resource Management**
   - Monitor memory usage
   - Adjust thread counts
   - Use appropriate containers

3. **Analysis Considerations**
   - Consider protein lengths
   - Check annotation quality
   - Validate clustering parameters

## 📚 Additional Resources

1. [Flye Documentation](https://github.com/fenderglass/Flye)
2. [Bakta Documentation](https://github.com/oschwengers/bakta)
3. [MMseqs2 User Guide](https://github.com/soedinglab/MMseqs2/wiki)
4. [Nextflow Patterns](https://github.com/nextflow-io/patterns)

## 🏆 Evaluation Criteria

1. Pipeline functionality (40%)
2. Code organization (20%)
3. Documentation quality (20%)
4. Analysis depth (20%)

Good luck with your project! Remember to commit your changes regularly and document your progress. 