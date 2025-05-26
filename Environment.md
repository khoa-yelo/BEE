# 🛠️ Setting Up Your Bioinformatics Environment

Welcome to the environment setup guide! In this section, we'll set up all the tools you need for bioinformatics work. We'll cover package management with Micromamba, containerization, and essential development tools.

## 📋 Table of Contents
- [Package Management with Micromamba](#package-management-with-micromamba)
- [Working with Containers](#working-with-containers)
- [Development Environments (IDEs)](#development-environments-ides)

## 📦 Package Management with Micromamba

### Understanding Package Management
Package management in Python scientific computing typically involves:
1. **Conda**: The traditional package manager for scientific computing
2. **Mamba**: A fast, C++-based reimplementation of Conda
3. **Micromamba**: A minimal, self-contained version of Mamba

### Why Micromamba?
- Faster than Conda and Mamba
- No base environment overhead
- Single binary, no dependencies
- Perfect for HPC environments

### Installation and Configuration
1. Install Micromamba:
   ```bash
   # Download and install
   curl -Ls https://micro.mamba.pm/api/micromamba/linux-64/latest | tar -xvj bin/micromamba
   
   # Move to a storage location with ample space (replace with your preferred path)
   mkdir -p /scratch/username/micromamba
   export MAMBA_ROOT_PREFIX=/scratch/username/micromamba
   
   # Initialize shell
   ./micromamba shell init -s bash -p $MAMBA_ROOT_PREFIX
   ```

2. Configure shell:
   ```bash
   # Add to your ~/.bashrc
   export MAMBA_ROOT_PREFIX=/scratch/username/micromamba
   eval "$(micromamba shell hook --shell bash)"
   ```

### Creating an Environment
```bash
# Create a new environment
micromamba create -n bioinfo python=3.10

# Activate the environment
micromamba activate bioinfo

# Install common bioinformatics packages
micromamba install -c bioconda -c conda-forge \
    biopython \
    pandas \
    numpy \
    matplotlib \
    seaborn \
    scikit-learn
```

## 🐳 Working with Containers

### Installing Docker Desktop

1. Download and install Docker Desktop:
   - Windows: [Download](https://desktop.docker.com/win/main/amd64/Docker%20Desktop%20Installer.exe)
   - Mac Intel: [Download](https://desktop.docker.com/mac/main/amd64/Docker.dmg)
   - Mac M1/M2: [Download](https://desktop.docker.com/mac/main/arm64/Docker.dmg)

2. Verify installation:
   ```bash
   docker --version
   docker run hello-world
   ```

### Setting Up Container Registry

1. Docker Hub (Free Public Registry):
   ```bash
   # Create account at https://hub.docker.com
   
   # Login via terminal
   docker login
   
   # Create repository on Docker Hub:
   # https://hub.docker.com/repository/create
   # Name: rnaseq
   # Visibility: Public/Private
   ```

2. GitHub Container Registry:
   ```bash
   # Create Personal Access Token (PAT)
   # Go to: GitHub → Settings → Developer settings → Personal access tokens
   # Select: write:packages, read:packages, delete:packages
   
   # Login to GitHub Container Registry
   export CR_PAT=YOUR_TOKEN
   echo $CR_PAT | docker login ghcr.io -u USERNAME --password-stdin
   
   # Tag image for GitHub
   docker tag yourusername/rnaseq:v1.0 ghcr.io/yourusername/rnaseq:v1.0
   ```

### Creating Docker Images

1. Create a Dockerfile:
   ```dockerfile
   # Base image with Python and R
   FROM mambaorg/micromamba:1.5.1
   
   # Set up working directory
   WORKDIR /app
   
   # Install packages using micromamba
   COPY environment.yml .
   RUN micromamba install -y -n base -f environment.yml && \
       micromamba clean -a -y
   
   # Install additional R packages for RNA-seq
   RUN R -e 'BiocManager::install(c( \
       "DESeq2", \
       "edgeR", \
       "tximport", \
       "ggplot2", \
       "pheatmap", \
       "EnhancedVolcano" \
   ), update=FALSE)'
   
   # Copy analysis scripts and data
   COPY scripts/ ./scripts/
   COPY data/ ./data/
   
   # Create output directory
   RUN mkdir output
   
   # Set default command (can be overridden)
   CMD ["Rscript", "scripts/rnaseq_analysis.R"]
   ```

2. Example `environment.yml`:
   ```yaml
   name: base
   channels:
     - conda-forge
     - bioconda
     - r
   dependencies:
     # Python packages
     - python=3.10
     - pandas
     - numpy
     - matplotlib
     - seaborn
     # R and Bioconductor
     - r-base=4.3
     - r-tidyverse
     - r-devtools
     - bioconductor-biocmanager
     - bioconductor-deseq2
     - bioconductor-edger
     - bioconductor-tximport
     # RNA-seq tools
     - salmon
     - fastqc
     - multiqc
     - star
   ```

3. Build and Push:
   ```bash
   # Build image
   docker build -t yourusername/rnaseq:v1.0 .
   
   # Login and push to registry
   docker login
   docker push yourusername/rnaseq:v1.0
   ```

### Using RNA-seq Container on HPC with Singularity

```bash
# Pull Docker image
singularity pull docker://yourusername/rnaseq:v1.0

# Run RNA-seq analysis
singularity exec \
  --bind ./data:/app/data \
  --bind ./output:/app/output \
  rnaseq_v1.0.sif Rscript scripts/rnaseq_analysis.R
```

## 💻 Development Environments (IDEs)

### Local Development Setup

1. VSCode Remote Connection:
   ```bash
   # 1. Install VSCode locally
   # Download from https://code.visualstudio.com/
   
   # 2. Install Remote Development extension pack
   
   # 3. Configure SSH (~/.ssh/config)
   Host hpc-server
       HostName your.hpc.address
       User your-username
       IdentityFile ~/.ssh/id_rsa
   
   # 4. Connect to Remote:
   # - Press F1 or Ctrl+Shift+P
   # - Type "Remote-SSH: Connect to Host"
   # - Select your configured host
   # - Wait for VSCode Server to install on remote
   
   # 5. Open folder on remote
   # - Click "Open Folder"
   # - Navigate to your project directory
   ```

### Server-Side Setup

1. Install packages:
   ```bash
   micromamba install -c conda-forge code-server jupyterlab nodejs
   ```

2. Start services:
   ```bash
   screen -S ide               # Create session
   code-server                 # Start VSCode
   Ctrl+a c                    # New window
   jupyter lab --no-browser   # Start Jupyter
   Ctrl+a d                    # Detach
   screen -r ide              # Return later
   ```

3. Access services:
   ```bash
   ssh -L 8080:localhost:8080 username@remote-server   # VSCode
   ssh -L 8888:localhost:8888 username@remote-server   # Jupyter
   ```

### Environment Management Tips

1. Storage Management:
   ```bash
   # Check environment sizes
   du -sh $MAMBA_ROOT_PREFIX/envs/*
   
   # Clean package cache
   micromamba clean --all
   ```

2. Extension Management:
   ```bash
   # For code-server, install extensions from CLI
   code-server --install-extension ms-python.python
   code-server --install-extension ms-toolsai.jupyter
   ```

3. Common Issues:
   - Port conflicts: Change ports in config files
   - Connection drops: Use tmux for session persistence
   - Memory issues: Monitor environment sizes and clean regularly
   - Extension loading: Check nodejs installation for VSCode extensions

## 🔍 Additional Resources

- [Micromamba Documentation](https://mamba.readthedocs.io/en/latest/user_guide/micromamba.html)
- [VSCode Remote Development](https://code.visualstudio.com/docs/remote/remote-overview)
- [Jupyter Documentation](https://jupyter.org/documentation)

