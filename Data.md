# 📊 Working with Biological Data

## Overview of Biological Data Types

### 🧬 Genomics
1. **DNA Sequencing**
   - Whole Genome Sequencing (WGS)
   - Exome Sequencing
   - ChIP-seq (protein-DNA interactions)
   - ATAC-seq (chromatin accessibility)
   - Hi-C (chromatin conformation)

2. **Variants**
   - Single Nucleotide Polymorphisms (SNPs)
   - Structural Variants
   - Copy Number Variations (CNVs)
   - Insertions/Deletions (Indels)

### 📝 Transcriptomics
1. **RNA Sequencing**
   - Bulk RNA-seq
   - Single-cell RNA-seq
   - Long-read RNA-seq
   - Small RNA sequencing
   - Ribosome profiling

2. **Expression Analysis**
   - Microarray data
   - Differential expression
   - Alternative splicing
   - Isoform quantification

### 🔬 Proteomics
1. **Protein Analysis**
   - Mass Spectrometry
   - Protein-protein interactions
   - Post-translational modifications
   - Protein structures (X-ray, Cryo-EM)

2. **Functional Analysis**
   - Protein activity assays
   - Protein localization
   - Pathway analysis

### 🧪 Metabolomics
1. **Small Molecules**
   - Metabolite profiles
   - Lipid profiles
   - Secondary metabolites

2. **Pathway Analysis**
   - Flux analysis
   - Metabolic networks
   - Biochemical pathways

### 🦠 Microbiome
1. **Community Analysis**
   - 16S rRNA sequencing
   - Metagenomic sequencing
   - Metatranscriptomics
   - Metabolomics

2. **Interaction Data**
   - Host-microbe interactions
   - Microbe-microbe interactions
   - Environmental factors

### 📸 Imaging
1. **Microscopy**
   - Light microscopy
   - Electron microscopy
   - Fluorescence microscopy
   - Live cell imaging
   - High-content screening

2. **Medical Imaging**
   - MRI
   - CT scans
   - X-rays
   - Ultrasound

### 📈 Clinical Data
1. **Patient Information**
   - Electronic Health Records
   - Clinical trials
   - Patient outcomes
   - Treatment responses

2. **Biomarkers**
   - Disease markers
   - Drug responses
   - Diagnostic indicators

### 🧮 Meta-Analysis
1. **Integration**
   - Multi-omics data
   - Cross-study comparisons
   - Meta-genomics

2. **Annotations**
   - Gene ontology
   - Pathway databases
   - Disease associations

## Fundamental Data Abstractions in Biology

When working with biological data programmatically, we can abstract most data types into these core computational representations:

### 1. Sequential Data (Ordered Series)
- **Definition**: Linear, ordered data where position matters
- **Examples**:
  - DNA/RNA sequences (ATCG strings)
  - Protein sequences (amino acid chains)
  - Time series (gene expression over time)
  - Evolutionary sequences/phylogenetic trees
- **Common Operations**:
  - Pattern matching
  - Alignment
  - Substring search
  - Sequential pattern mining
- **Typical Formats**: FASTA, FASTQ, Stockholm

### 2. Matrix/Tensor Data (N-dimensional Arrays)
- **Definition**: Data organized in multiple dimensions
- **Examples**:
  - Gene expression matrices (genes × samples)
  - Single-cell data (cells × genes)
  - Imaging data (x × y × z × channels)
  - Contact matrices (Hi-C data)
- **Common Operations**:
  - Dimensionality reduction
  - Clustering
  - Statistical analysis
  - Image processing
- **Typical Formats**: HDF5, NetCDF, TIFF

### 3. Tabular Data (Relational)
- **Definition**: Data organized in rows and columns with relationships
- **Examples**:
  - Sample metadata
  - Clinical data
  - Variant annotations
  - Pathway databases
- **Common Operations**:
  - Filtering
  - Joining
  - Aggregation
  - Statistical summaries
- **Typical Formats**: CSV, TSV, SQL databases

### 4. Graph/Network Data
- **Definition**: Data representing connections between entities
- **Examples**:
  - Protein-protein interaction networks
  - Metabolic pathways
  - Gene regulatory networks
  - Phylogenetic trees
- **Common Operations**:
  - Path finding
  - Network analysis
  - Community detection
  - Centrality measures
- **Typical Formats**: GraphML, JSON, Neo4j

### 5. Key-Value/Document Data
- **Definition**: Schemaless, flexible data structures
- **Examples**:
  - Experimental metadata
  - Annotations
  - Configuration data
  - Unstructured clinical notes
- **Common Operations**:
  - Flexible querying
  - Text mining
  - Document classification
  - Information extraction
- **Typical Formats**: JSON, YAML, MongoDB

### Cross-cutting Considerations
1. **Scale**
   - Small (fits in memory)
   - Medium (needs disk)
   - Big (distributed processing)

2. **Access Patterns**
   - Random access vs. Sequential
   - Read-heavy vs. Write-heavy
   - Batch vs. Streaming

3. **Quality Control**
   - Data validation
   - Error checking
   - Version control
   - Provenance tracking

## Practical Implementations of Data Abstractions

### 1. Working with Sequential Data
```python
import sqlite3
import Bio.SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

class SequenceHandler:
    def __init__(self, db_path):
        self.db_path = db_path
        self._init_db()
    
    def _init_db(self):
        """Initialize SQLite database with sequence-specific indices"""
        with sqlite3.connect(self.db_path) as conn:
            conn.execute('''
                CREATE TABLE IF NOT EXISTS sequences (
                    seq_id TEXT PRIMARY KEY,
                    sequence TEXT,
                    length INTEGER,
                    metadata JSON,
                    UNIQUE(seq_id)
                )
            ''')
            # Index for sequence pattern matching
            conn.execute('CREATE INDEX IF NOT EXISTS seq_pattern ON sequences(sequence)')
    
    def store_sequence(self, seq_record):
        """Store a sequence with its metadata"""
        with sqlite3.connect(self.db_path) as conn:
            conn.execute(
                'INSERT INTO sequences VALUES (?, ?, ?, ?)',
                (seq_record.id, 
                 str(seq_record.seq),
                 len(seq_record.seq),
                 json.dumps(seq_record.annotations))
            )
    
    def find_pattern(self, pattern):
        """Find sequences containing a specific pattern"""
        with sqlite3.connect(self.db_path) as conn:
            return conn.execute(
                'SELECT seq_id, sequence FROM sequences WHERE sequence LIKE ?',
                (f'%{pattern}%',)
            ).fetchall()

# Example usage
handler = SequenceHandler('sequences.db')
seq_record = SeqRecord(
    Seq("ATCGATCG"),
    id="seq1",
    annotations={"organism": "Human", "gene": "BRCA1"}
)
handler.store_sequence(seq_record)
```

### 2. Working with Matrix Data
```python
import h5py
import numpy as np
from typing import Tuple, Optional

class MatrixDataset:
    def __init__(self, filename: str):
        self.filename = filename
        
    def create_dataset(self, 
                      name: str, 
                      data: np.ndarray,
                      chunk_size: Optional[Tuple] = None):
        """Store matrix data with automatic chunking"""
        with h5py.File(self.filename, 'a') as f:
            if chunk_size is None:
                # Automatically determine chunk size
                chunk_size = self._optimal_chunk_size(data.shape)
            
            # Create dataset with compression
            f.create_dataset(
                name,
                data=data,
                chunks=chunk_size,
                compression='gzip',
                compression_opts=4
            )
            
            # Store metadata
            f[name].attrs['shape'] = data.shape
            f[name].attrs['created'] = np.string_(datetime.now().isoformat())
    
    def _optimal_chunk_size(self, shape: Tuple) -> Tuple:
        """Calculate optimal chunk size based on data shape"""
        chunk_size = list(shape)
        for i in range(len(chunk_size)):
            chunk_size[i] = min(chunk_size[i], 1000)
        return tuple(chunk_size)
    
    def get_slice(self, name: str, slices: Tuple):
        """Efficiently retrieve a slice of data"""
        with h5py.File(self.filename, 'r') as f:
            return f[name][slices]

# Example usage
matrix_handler = MatrixDataset('expression_data.h5')
expression_matrix = np.random.rand(10000, 1000)  # genes x samples
matrix_handler.create_dataset('expression', expression_matrix)
```

### 3. Working with Relational Data
```python
import sqlite3
import json
from typing import Dict, Optional, List
from datetime import datetime

class BiologicalDatabase:
    def __init__(self, db_path: str):
        self.db_path = db_path
        self._init_db()
    
    def _init_db(self):
        """Initialize the database schema"""
        with sqlite3.connect(self.db_path) as conn:
            # Samples table
            conn.execute('''
                CREATE TABLE IF NOT EXISTS samples (
                    sample_id TEXT PRIMARY KEY,
                    condition TEXT NOT NULL,
                    collection_date TEXT,
                    created_at TEXT DEFAULT CURRENT_TIMESTAMP
                )
            ''')
            
            # Sample metadata table (key-value pairs)
            conn.execute('''
                CREATE TABLE IF NOT EXISTS sample_metadata (
                    sample_id TEXT,
                    key TEXT,
                    value TEXT,
                    PRIMARY KEY (sample_id, key),
                    FOREIGN KEY (sample_id) REFERENCES samples(sample_id)
                )
            ''')
            
            # Measurements table (e.g., gene expression, protein levels)
            conn.execute('''
                CREATE TABLE IF NOT EXISTS measurements (
                    id INTEGER PRIMARY KEY AUTOINCREMENT,
                    sample_id TEXT,
                    feature_id TEXT,
                    value REAL,
                    measurement_type TEXT,
                    FOREIGN KEY (sample_id) REFERENCES samples(sample_id)
                )
            ''')
            
            # Create indices for common queries
            conn.execute('CREATE INDEX IF NOT EXISTS idx_sample_condition ON samples(condition)')
            conn.execute('CREATE INDEX IF NOT EXISTS idx_measurement_feature ON measurements(feature_id)')
    
    def add_sample(self, sample_id: str, condition: str, 
                  metadata: Optional[Dict] = None,
                  collection_date: Optional[str] = None):
        """Add a new sample with optional metadata"""
        with sqlite3.connect(self.db_path) as conn:
            try:
                # Insert sample
                conn.execute(
                    'INSERT INTO samples (sample_id, condition, collection_date) VALUES (?, ?, ?)',
                    (sample_id, condition, collection_date)
                )
                
                # Insert metadata if provided
                if metadata:
                    metadata_records = [
                        (sample_id, key, str(value))
                        for key, value in metadata.items()
                    ]
                    conn.executemany(
                        'INSERT INTO sample_metadata (sample_id, key, value) VALUES (?, ?, ?)',
                        metadata_records
                    )
                
            except sqlite3.IntegrityError as e:
                raise ValueError(f"Sample {sample_id} already exists or invalid data: {e}")
    
    def add_measurements(self, sample_id: str, measurements: List[Dict]):
        """Add measurements for a sample"""
        with sqlite3.connect(self.db_path) as conn:
            try:
                records = [
                    (sample_id, m['feature_id'], m['value'], m.get('type', 'unknown'))
                    for m in measurements
                ]
                conn.executemany(
                    '''INSERT INTO measurements 
                       (sample_id, feature_id, value, measurement_type)
                       VALUES (?, ?, ?, ?)''',
                    records
                )
            except sqlite3.IntegrityError as e:
                raise ValueError(f"Invalid measurement data for sample {sample_id}: {e}")
    
    def get_sample_data(self, sample_id: str) -> Dict:
        """Get all data for a sample including metadata and measurements"""
        with sqlite3.connect(self.db_path) as conn:
            # Get sample info
            sample = conn.execute(
                'SELECT * FROM samples WHERE sample_id = ?',
                (sample_id,)
            ).fetchone()
            
            if not sample:
                raise ValueError(f"Sample {sample_id} not found")
            
            # Get metadata
            metadata = dict(conn.execute(
                'SELECT key, value FROM sample_metadata WHERE sample_id = ?',
                (sample_id,)
            ).fetchall())
            
            # Get measurements
            measurements = conn.execute(
                'SELECT feature_id, value, measurement_type FROM measurements WHERE sample_id = ?',
                (sample_id,)
            ).fetchall()
            
            return {
                'sample_id': sample[0],
                'condition': sample[1],
                'collection_date': sample[2],
                'created_at': sample[3],
                'metadata': metadata,
                'measurements': [
                    {'feature_id': m[0], 'value': m[1], 'type': m[2]}
                    for m in measurements
                ]
            }
    
    def query_samples(self, condition: Optional[str] = None) -> List[Dict]:
        """Query samples with optional condition filter"""
        with sqlite3.connect(self.db_path) as conn:
            query = 'SELECT sample_id, condition, collection_date FROM samples'
            params = []
            
            if condition:
                query += ' WHERE condition = ?'
                params.append(condition)
            
            return [
                {
                    'sample_id': row[0],
                    'condition': row[1],
                    'collection_date': row[2]
                }
                for row in conn.execute(query, params).fetchall()
            ]
    
    def get_sample_statistics(self, sample_id: str) -> Dict:
        """Get statistical summary of sample measurements"""
        with sqlite3.connect(self.db_path) as conn:
            return conn.execute('''
                SELECT 
                    COUNT(*) as count,
                    AVG(value) as mean,
                    MIN(value) as min,
                    MAX(value) as max,
                    measurement_type
                FROM measurements 
                WHERE sample_id = ?
                GROUP BY measurement_type
            ''', (sample_id,)).fetchall()

# Example usage
db = BiologicalDatabase('bio_samples.db')

# Add a sample with metadata
db.add_sample(
    sample_id='S1',
    condition='control',
    metadata={
        'age': 25,
        'tissue': 'liver',
        'treatment': 'none'
    },
    collection_date='2024-03-15'
)

# Add measurements
db.add_measurements('S1', [
    {'feature_id': 'GENE1', 'value': 123.45, 'type': 'expression'},
    {'feature_id': 'GENE2', 'value': 67.89, 'type': 'expression'},
    {'feature_id': 'PROTEIN1', 'value': 0.23, 'type': 'abundance'}
])

# Query data
sample_data = db.get_sample_data('S1')
control_samples = db.query_samples(condition='control')
statistics = db.get_sample_statistics('S1')
```

### 4. Working with Graph Data
```python
import networkx as nx
import json

class BiologicalNetwork:
    def __init__(self):
        self.graph = nx.Graph()
    
    def add_interaction(self, source, target, interaction_type, score):
        """Add a biological interaction to the network"""
        self.graph.add_edge(
            source, 
            target, 
            type=interaction_type,
            score=score
        )
    
    def find_modules(self, min_size=3):
        """Identify biological modules/communities"""
        communities = nx.community.louvain_communities(self.graph)
        return [c for c in communities if len(c) >= min_size]
    
    def export_graph(self, filename):
        """Export network in standard format"""
        data = nx.node_link_data(self.graph)
        with open(filename, 'w') as f:
            json.dump(data, f)

# Example usage
network = BiologicalNetwork()
network.add_interaction('BRCA1', 'BRCA2', 'protein_interaction', 0.9)
network.add_interaction('BRCA1', 'TP53', 'genetic_interaction', 0.7)
```

### 5. Working with Document Data
```python
from pymongo import MongoClient
from datetime import datetime

class ExperimentalMetadata:
    def __init__(self, connection_string):
        self.client = MongoClient(connection_string)
        self.db = self.client.biodata
    
    def store_experiment(self, experiment_data):
        """Store flexible experimental metadata"""
        experiment_data['timestamp'] = datetime.now()
        return self.db.experiments.insert_one(experiment_data)
    
    def find_experiments(self, query):
        """Flexible querying of experimental data"""
        return self.db.experiments.find(query)
    
    def update_experiment(self, experiment_id, updates):
        """Update experimental metadata"""
        return self.db.experiments.update_one(
            {'_id': experiment_id},
            {'$set': updates}
        )

# Example usage
metadata_handler = ExperimentalMetadata('mongodb://localhost:27017/')
experiment = {
    'name': 'RNA-seq_001',
    'protocol': {'type': 'single-cell', 'version': '2.0'},
    'conditions': ['control', 'treatment'],
    'notes': 'Sample quality looks good'
}
metadata_handler.store_experiment(experiment)
```

## Best Practices for Data Integration

### 1. Data Validation
```python
from pydantic import BaseModel, Field
from typing import List, Dict, Optional

class ExperimentMetadata(BaseModel):
    experiment_id: str
    description: str
    date: datetime
    conditions: List[str]
    replicates: int
    notes: Optional[str] = None

# Validate data before storage
metadata = ExperimentMetadata(
    experiment_id="EXP001",
    description="RNA-seq analysis",
    date=datetime.now(),
    conditions=["control", "treatment"],
    replicates=3
)
```

### 2. Data Integration Example
```python
def integrate_multi_omics(sequence_db, expression_h5, metadata_db, sample_id):
    """Integrate data across different storage systems"""
    # Get sequence data
    seq_handler = SequenceHandler(sequence_db)
    sequences = seq_handler.get_sequences(sample_id)
    
    # Get expression data
    matrix_handler = MatrixDataset(expression_h5)
    expression = matrix_handler.get_slice('expression', sample_id)
    
    # Get metadata
    bio_db = BiologicalDatabase(metadata_db)
    metadata = bio_db.get_sample_data(sample_id)
    
    return {
        'sequences': sequences,
        'expression': expression,
        'metadata': metadata
    }
```

Remember: Choose the right data abstraction and storage solution for your specific use case! 🎯 