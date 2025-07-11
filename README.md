# LAMPforge - Advanced LAMP Primer Designer

[![Python 3.8+](https://img.shields.io/badge/python-3.8+-blue.svg)](https://www.python.org/downloads/)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

LAMPforge is a comprehensive Python tool for designing Loop-mediated isothermal amplification (LAMP) primers with advanced features including machine learning scoring, type-specific thresholds, and batch processing capabilities.

## 🧬 Features

- **Flexible Primer Design**: Supports various sequence types (normal, bisulfite, pyrosequencing)
- **Progressive Constraint Relaxation**: Automatically adapts to challenging sequences (AT-rich, GC-rich)
- **Type-Specific Thresholds**: Optimizes parameters for each primer type (F3, B3, FIP, BIP, LF, LB)
- **Machine Learning Scoring**: Optional ML-based primer set ranking
- **Batch Processing**: Design primers for multiple targets simultaneously
- **Loop Primer Design**: Automated LF/LB primer design with proper topology understanding
- **Specificity Checking**: K-mer based off-target screening
- **Multiple Output Formats**: CSV, JSON, and detailed analysis reports

## 📋 Requirements

```bash
python >= 3.8
numpy >= 1.19.0
pandas >= 1.3.0
```

Optional dependencies:
```bash
scikit-learn >= 0.24.0  # For ML scoring
scipy >= 1.7.0 # For scientific computing and math operations
biopython >= 1.79 # For advanced biology data processing
tqdm >= 4.62.0  # Progress bars
```

Visualization dependencies:
```bash
matplotlib >= 3.4.0
seaborn >= 0.11.0
```

## 🚀 Installation

### Clone the repository
```bash
git clone https://github.com/yourusername/lampforge.git
cd lampforge
```

### Install dependencies
```bash
pip install -r requirements.txt
```

### Or install with conda
```bash
conda env create -f environment.yml
conda activate lampforge
```

## 💻 Usage

### Basic Usage

```bash
python lamp_designer.py -i sequence.fasta -o results/
```

### Direct sequence input
```bash
python lamp_designer.py -i "ATGCGATCGTAGCTAGCTAGCGATCGTAGCTAGC..." -o results/
```

### Advanced Options

```bash
# Design with loop primers
python lamp_designer.py -i sequence.fasta -o results/ --loop_primers

# Use machine learning scoring
python lamp_designer.py -i sequence.fasta -o results/ --ml_scoring

# Disable quality filtering (get ALL combinations)
python lamp_designer.py -i sequence.fasta -o results/ --full

# Custom parameters
python lamp_designer.py -i sequence.fasta -o results/ \
    -t 55 -T 75 \          # Tm range: 55-75°C
    -c 30 -C 70 \          # GC content: 30-70%
    -n 5                   # Return top 5 primer sets
```

### Type-Specific Design (Experimental)

```bash
# Enable primer-type-specific threshold optimization
python lamp_designer.py -i sequence.fasta -o results/ \
    --enable-type-specific \
    --print-type-analysis
```

## 📊 Input Formats

### FASTA file
```
>Sequence_name
ATGCGATCGTAGCTAGCTAGCGATCGTAGCTAGCGATCGTAGCTAGCTAGCGATCG
TAGCTAGCGATCGTAGCTAGCTAGCGATCGTAGCTAGCGATCGTAGCTAGCTAGCG
```

### Direct sequence
```bash
python lamp_designer.py -i "ATGCGATCGTAGCTAGCTAGCGATCG" -o results/
```

## 📁 Output Files

The tool generates multiple output files in the specified directory:

- `lamp_primers_[timestamp]_top.csv` - Top primer sets in CSV format
- `lamp_primers_[timestamp]_top.json` - Top primer sets in JSON format
- `lamp_primers_[timestamp]_analysis.txt` - Detailed analysis report
- `lamp_primers_[timestamp]_config.json` - Design parameters used
- `synthesis_order.json` - Ready-to-order primer sequences (if type-specific design)

### Example Output (CSV)

```csv
Set_ID,F3_Sequence,F3_Tm,F3_GC,B3_Sequence,B3_Tm,B3_GC,FIP_Sequence,FIP_Tm,FIP_GC,BIP_Sequence,BIP_Tm,BIP_GC
Set_001,GCGATCGTAGCTAGCTA,62.5,52.9,TAGCGATCGTAGCTAG,61.8,50.0,CTAGCTAGCGATCGTAGCTAGCGATCGTAGCTA,63.2,51.5,GCTAGCGATCGTAGCTAGCTAGCGATCGTAG,62.9,54.8
```

## 🔧 Command Line Options

### Required Arguments
- `-i, --input`: Input sequence (FASTA file or string)
- `-o, --output`: Output directory for results

### Optional Arguments

#### Output Options
- `-f, --format`: Result format: `top_n` or `detailed` (default: top_n)
- `-n, --top_n`: Number of top primer sets to return (default: 3)
- `--full`: Save ALL primer combinations (disables quality filtering)

#### Primer Constraints
- `-l, --min_length`: Minimum primer length (default: 15)
- `-L, --max_length`: Maximum primer length (default: 40)
- `-t, --min_tm`: Minimum melting temperature (default: 50.0°C)
- `-T, --max_tm`: Maximum melting temperature (default: 80.0°C)
- `-c, --min_gc`: Minimum GC content % (default: 20.0)
- `-C, --max_gc`: Maximum GC content % (default: 80.0)

#### Thermodynamic Constraints
- `-g, --max_3prime_dg`: Maximum ΔG for 3' end stability (default: -0.5 kcal/mol)
- `-G, --max_hairpin_dg`: Maximum ΔG for hairpin formation (default: -1.0 kcal/mol)
- `-d, --max_dimer_dg`: Maximum ΔG for primer-dimer formation (default: -3.0 kcal/mol)

#### Advanced Options
- `--ml_scoring`: Enable machine learning scoring
- `--loop_primers`: Design loop primers (LF/LB)
- `--no_specificity`: Disable specificity checking
- `-s, --seq_type`: Sequence type: `normal`, `bisulfite`, or `pyrosequencing`

## 🧪 Examples

### Example 1: Basic LAMP primer design
```bash
python lamp_designer.py -i examples/target1.fasta -o results/target1/
```

### Example 2: AT-rich sequence with relaxed constraints
```bash
python lamp_designer.py -i examples/at_rich.fasta -o results/at_rich/ \
    -t 45 -T 75 \      # Wider Tm range
    -c 15 -C 60        # Lower GC content allowed
```

### Example 3: Full analysis with all features
```bash
python lamp_designer.py -i examples/target2.fasta -o results/full_analysis/ \
    --loop_primers \
    --ml_scoring \
    -f detailed \
    -n 10
```

## 🔬 Algorithm Details

PyLAMP uses a multi-stage approach for primer design:

1. **Sequence Analysis**: Analyzes GC content and composition
2. **Progressive Search**: Applies increasingly relaxed constraints to find candidates
3. **Type-Specific Optimization**: Adjusts parameters based on primer type
4. **Tm-Balanced Selection**: Ensures uniform melting temperatures across primer sets
5. **Quality Scoring**: Ranks primer sets based on multiple criteria
6. **Loop Primer Addition**: Designs LF/LB primers for enhanced amplification

### Scoring Criteria

- Tm uniformity across all primers
- GC content balance
- 3' end stability
- Minimal secondary structure formation
- Proper primer positioning for LAMP topology

## 📈 Performance Tips

1. **For difficult sequences**: Use `--full` mode to see all possible combinations
2. **For AT-rich sequences**: Lower minimum Tm and GC constraints
3. **For GC-rich sequences**: Increase maximum Tm and adjust GC range
4. **For short sequences**: Reduce minimum primer lengths

## 🐛 Troubleshooting

### No primers found
- Try relaxing constraints (wider Tm range, GC content)
- Use `--full` mode to disable filtering
- Check sequence length (minimum ~200bp recommended)

### Too many primers
- Tighten constraints
- Reduce `-n` parameter
- Enable specificity checking

## 📚 Citation

If you use LAMPforge in your research, please cite:

```
LAMPforge: Advanced LAMP Primer Designer
Sukhanova Xenia, 2025
https://github.com/sukhanovaxenia/lampforge
```

## 🤝 Contributing

Contributions are welcome! Please feel free to submit a Pull Request. For major changes, please open an issue first to discuss what you would like to change.

## 📄 License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.

## 🙏 Acknowledgments

- LAMP primer design principles based on Notomi et al. (2000)
- Thermodynamic calculations adapted from SantaLucia (1998)
- Nearest neighbor parameters from Allawi & SantaLucia (1997)

## 📮 Contact

For questions, issues, or suggestions, please open an issue on GitHub or contact sukhanovaxenia@gmail.com

---

**Note**: This tool is for research use only. Always validate primer designs experimentally before use.
