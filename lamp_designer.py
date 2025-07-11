#!/usr/bin/env python3
"""
LAMP Primer Designer - Core Module
A research tool for designing LAMP (Loop-mediated isothermal amplification) primers
"""

import re
import os
import json
import csv
import math
from argparse import ArgumentParser
from typing import Dict, List, Tuple, Optional, NamedTuple, Union, Set
from dataclasses import dataclass, asdict
from enum import Enum
import itertools
import numpy as np
import hashlib
from collections import defaultdict
from datetime import datetime

class SequenceType(Enum):
    NORMAL = "normal"
    BISULFITE = "bisulfite"
    PYROSEQUENCING = "pyrosequencing"

class ResultFormat(Enum):
    TOP_N = "top_n"
    ALL = "all"
    DETAILED = "detailed"

@dataclass
class PrimerConfig:
    """Configuration for LAMP primer design"""
    # Primer length constraints
    min_length: int = 15          # Reduced from 18
    max_length: int = 40          # Increased from 25
    
    # Melting temperature constraints (°C)
    min_tm: float = 30.0          # Reduced from 58.0
    max_tm: float = 90.0          # Increased from 65.0
    
    # GC content constraints (%) - CRITICAL FIX
    min_gc: float = 10.0          # Reduced from 40.0 (was blocking low-GC primers)
    max_gc: float = 90.0          # Added upper limit (was missing)
    
    # Distance constraints (bp) - kept same
    min_f1c_f2c_distance: int = 10
    max_f1c_f2c_distance: int = 80
    min_b1c_b2c_distance: int = 10
    max_b1c_b2c_distance: int = 80
    
    # Loop primer constraints
    loop_primer_length: int = 20
    loop_distance_from_inner: int = 20
    
    # Gibbs free energy constraints (kcal/mol) - More permissive
    max_3prime_dg: float = -0.5   # Less strict (was -2.0)
    max_hairpin_dg: float = -1.0  # Less strict (was -2.0)
    max_dimer_dg: float = -3.0    # Less strict (was -5.0)
    
    # Reaction conditions - kept same
    na_conc: float = 50.0  # mM
    k_conc: float = 50.0   # mM
    mg_conc: float = 8.0   # mM
    dntp_conc: float = 0.8 # mM
    temperature: float = 65.0  # Reaction temperature (°C)
    
    # Non-specificity control
    enable_specificity_check: bool = True
    max_off_targets: int = 3
    min_off_target_tm: float = 45.0
    specificity_database: Optional[str] = None
    kmer_size: int = 15
    
    # Sequence type
    seq_type: SequenceType = SequenceType.NORMAL
    
    # Use ML for scoring if enabled
    use_ml_scoring: bool = False
    ml_training_size: int = 1000

class Primer(NamedTuple):
    """Represents a primer with its properties"""
    sequence: str
    start: int
    end: int
    tm: float
    gc_content: float
    primer_type: str  # F3, B3, FIP, BIP, LF, LB
    dg_3prime: float  # ΔG at 3' end (kcal/mol)
    dg_hairpin: float  # ΔG for hairpin formation (kcal/mol)

class LAMPPrimerSet(NamedTuple):
    """Complete set of LAMP primers"""
    f3: Primer
    b3: Primer
    fip: Primer  # F1c + F2
    bip: Primer  # B1c + B2
    lf: Optional[Primer] = None
    lb: Optional[Primer] = None

class LAMPPrimerType(Enum):
    """LAMP primer types with their functional roles"""
    F3 = "F3"      # Forward outer primer
    B3 = "B3"      # Backward outer primer  
    F1C = "F1c"    # Forward inner complement (part of FIP)
    F2 = "F2"      # Forward inner primer (part of FIP)
    B1C = "B1c"    # Backward inner complement (part of BIP)
    B2 = "B2"      # Backward inner primer (part of BIP)
    FIP = "FIP"    # Forward inner primer (composite F1c + F2)
    BIP = "BIP"    # Backward inner primer (composite B1c + B2)
    LF = "LF"      # Loop forward primer
    LB = "LB"      # Loop backward primer

@dataclass 
class PrimerTypeSpecificConfig:
    """Configuration with primer-type-specific thresholds"""
    
    # Primer-type-specific parameter ranges
    primer_type_configs: Dict[LAMPPrimerType, Dict] = None
    
    # Base composition adjustments  
    composition_adjustments: Dict = None
    
    def __post_init__(self):
        """Initialize primer-type-specific configurations"""
        if self.primer_type_configs is None:
            self.primer_type_configs = {
                
                # OUTER PRIMERS (F3/B3) - MUCH more permissive for real sequences
                LAMPPrimerType.F3: {
                    "tm_range": (40.0, 85.0),          # MUCH wider (was 58-72)
                    "tm_optimal": 65.0,
                    "length_range": (15, 30),          # More flexible
                    "gc_range": (10.0, 90.0),          # MUCH wider (was 30-70)
                    "dg_3prime_max": 5.0,              # VERY permissive (was -0.8)
                    "dg_hairpin_max": 5.0,             # VERY permissive (was -1.5)
                    "dg_dimer_max": 5.0,               # VERY permissive (was -4.0)
                    "priority_weight": 1.2,
                    "description": "Forward outer primer - PERMISSIVE"
                },
                
                LAMPPrimerType.B3: {
                    "tm_range": (40.0, 85.0),          
                    "tm_optimal": 65.0,
                    "length_range": (15, 30),          
                    "gc_range": (10.0, 90.0),          
                    "dg_3prime_max": 5.0,              
                    "dg_hairpin_max": 5.0,             
                    "dg_dimer_max": 5.0,               
                    "priority_weight": 1.2,        
                    "description": "Backward outer primer - PERMISSIVE"
                },
                
                # INNER PRIMER COMPONENTS (F1c/B1c) - VERY permissive
                LAMPPrimerType.F1C: {
                    "tm_range": (35.0, 85.0),          # VERY wide range
                    "tm_optimal": 62.0,
                    "length_range": (12, 25),          # Wider range
                    "gc_range": (5.0, 95.0),           # EXTREMELY wide GC range
                    "dg_3prime_max": 10.0,             # EXTREMELY permissive
                    "dg_hairpin_max": 10.0,            # EXTREMELY permissive
                    "dg_dimer_max": 10.0,              # EXTREMELY permissive
                    "priority_weight": 1.0,
                    "description": "F1 complement - VERY PERMISSIVE"
                },
                
                LAMPPrimerType.B1C: {
                    "tm_range": (35.0, 85.0),          # Same as F1c
                    "tm_optimal": 62.0,
                    "length_range": (12, 25),          
                    "gc_range": (5.0, 95.0),           
                    "dg_3prime_max": 10.0,             
                    "dg_hairpin_max": 10.0,            
                    "dg_dimer_max": 10.0,              
                    "priority_weight": 1.0,        
                    "description": "B1 complement - VERY PERMISSIVE"
                },
                
                # INNER PRIMER SEQUENCES (F2/B2) - VERY permissive
                LAMPPrimerType.F2: {
                    "tm_range": (35.0, 85.0),          # VERY wide range
                    "tm_optimal": 63.0,
                    "length_range": (12, 25),          
                    "gc_range": (5.0, 95.0),           
                    "dg_3prime_max": 10.0,             
                    "dg_hairpin_max": 10.0,            
                    "dg_dimer_max": 10.0,              
                    "priority_weight": 1.1,
                    "description": "Forward inner primer - VERY PERMISSIVE"
                },
                
                LAMPPrimerType.B2: {
                    "tm_range": (35.0, 85.0),          # Same as F2
                    "tm_optimal": 63.0,
                    "length_range": (12, 25),          
                    "gc_range": (5.0, 95.0),           
                    "dg_3prime_max": 10.0,             
                    "dg_hairpin_max": 10.0,            
                    "dg_dimer_max": 10.0,              
                    "priority_weight": 1.1,        
                    "description": "Backward inner primer - VERY PERMISSIVE"
                },
                
                # COMPOSITE PRIMERS (FIP/BIP) - VERY permissive
                LAMPPrimerType.FIP: {
                    "tm_range": (35.0, 85.0),          # VERY wide range
                    "tm_optimal": 62.5,
                    "length_range": (25, 60),          # Longer max length
                    "gc_range": (10.0, 90.0),          
                    "dg_3prime_max": 10.0,             
                    "dg_hairpin_max": 10.0,            
                    "dg_dimer_max": 10.0,              
                    "priority_weight": 1.5,
                    "description": "Forward inner composite - VERY PERMISSIVE"
                },
                
                LAMPPrimerType.BIP: {
                    "tm_range": (35.0, 85.0),          # Same as FIP
                    "tm_optimal": 62.5,
                    "length_range": (25, 60),          
                    "gc_range": (10.0, 90.0),          
                    "dg_3prime_max": 10.0,             
                    "dg_hairpin_max": 10.0,            
                    "dg_dimer_max": 10.0,              
                    "priority_weight": 1.5,        
                    "description": "Backward inner composite - VERY PERMISSIVE"
                },
                
                # LOOP PRIMERS (LF/LB) - EXTREMELY permissive
                LAMPPrimerType.LF: {
                    "tm_range": (30.0, 90.0),          # EXTREMELY wide range
                    "tm_optimal": 60.0,
                    "length_range": (12, 30),          
                    "gc_range": (0.0, 100.0),          # EXTREMELY wide GC range
                    "dg_3prime_max": 20.0,             # EXTREMELY permissive
                    "dg_hairpin_max": 20.0,            
                    "dg_dimer_max": 20.0,              
                    "priority_weight": 0.8,
                    "description": "Loop forward - EXTREMELY PERMISSIVE"
                },
                
                LAMPPrimerType.LB: {
                    "tm_range": (30.0, 90.0),          # Same as LF
                    "tm_optimal": 60.0,
                    "length_range": (12, 30),          
                    "gc_range": (0.0, 100.0),          
                    "dg_3prime_max": 20.0,             
                    "dg_hairpin_max": 20.0,            
                    "dg_dimer_max": 20.0,              
                    "priority_weight": 0.8,        
                    "description": "Loop backward - EXTREMELY PERMISSIVE"
                }
            }
    
        if self.composition_adjustments is None:
            # MUCH more permissive composition adjustments
            self.composition_adjustments = {
                "AT_RICH": {  # GC < 40%
                    "tm_adjustment": -1.0,     # Smaller adjustment (was -3.0)
                    "gc_tolerance": 30.0,      # MUCH more tolerance (was 15.0)
                    "length_adjustment": 1,     # Smaller adjustment (was 2)
                    "dg_relaxation": 3.0       # MUCH more permissive ΔG (was 1.5)
                },
                "NORMAL": {   # GC 40-65%
                    "tm_adjustment": 0.0,      # No adjustment
                    "gc_tolerance": 25.0,      # More tolerance (was 0.0)
                    "length_adjustment": 0,     # Standard length
                    "dg_relaxation": 3.0       # MUCH more permissive ΔG (was 1.0)
                },
                "GC_RICH": {  # GC > 65%
                    "tm_adjustment": 1.0,      # Smaller adjustment (was 2.0)
                    "gc_tolerance": 20.0,      # Less strict (was -5.0)
                    "length_adjustment": 0,     # No length penalty (was -1)
                    "dg_relaxation": 3.0       # MUCH more permissive ΔG (was 0.8)
                }
            }

class PrimerTypeSpecificDesigner:
    """LAMP designer with primer-type-specific parameter management"""
    
    def __init__(self, config: PrimerTypeSpecificConfig):
        self.config = config
        self.target_composition = None
        self.adapted_configs = {}
        
    def classify_target_composition(self, sequence: str) -> str:
        """Classify target as AT_RICH, NORMAL, or GC_RICH"""
        gc_content = (sequence.count('G') + sequence.count('C')) / len(sequence) * 100
        
        if gc_content < 40:
            return "AT_RICH"
        elif gc_content > 65:
            return "GC_RICH"
        else:
            return "NORMAL"
    
    def get_adapted_config_for_primer_type(self, primer_type: LAMPPrimerType, 
                                         composition: str) -> Dict:
        """Get composition-adapted configuration for specific primer type"""
        
        base_config = self.config.primer_type_configs[primer_type].copy()
        composition_adj = self.config.composition_adjustments[composition]
        
        # Apply composition adjustments
        tm_min, tm_max = base_config["tm_range"]
        adapted_config = {
            "tm_range": (
                tm_min + composition_adj["tm_adjustment"],
                tm_max + composition_adj["tm_adjustment"]
            ),
            "tm_optimal": base_config["tm_optimal"] + composition_adj["tm_adjustment"],
            
            "length_range": (
                max(10, base_config["length_range"][0] + composition_adj["length_adjustment"]),
                base_config["length_range"][1] + composition_adj["length_adjustment"]
            ),
            
            "gc_range": (
                max(5.0, base_config["gc_range"][0] - composition_adj["gc_tolerance"]),
                min(95.0, base_config["gc_range"][1] + composition_adj["gc_tolerance"])
            ),
            
            "dg_3prime_max": base_config["dg_3prime_max"] * composition_adj["dg_relaxation"],
            "dg_hairpin_max": base_config["dg_hairpin_max"] * composition_adj["dg_relaxation"],
            "dg_dimer_max": base_config["dg_dimer_max"] * composition_adj["dg_relaxation"],
            
            "priority_weight": base_config["priority_weight"],
            "description": base_config["description"],
            "composition": composition,
            "primer_type": primer_type.value
        }
        
        return adapted_config
    
    def analyze_and_adapt_all_primers(self, sequence: str) -> Dict[LAMPPrimerType, Dict]:
        """Analyze sequence and create adapted configs for all primer types"""
        
        self.target_composition = self.classify_target_composition(sequence)
        gc_content = (sequence.count('G') + sequence.count('C')) / len(sequence) * 100
        
        print(f"\n🧬 TARGET SEQUENCE ANALYSIS:")
        print(f"{'='*50}")
        print(f"📏 Length: {len(sequence)} bp")
        print(f"📊 GC content: {gc_content:.1f}%")
        print(f"🎯 Classification: {self.target_composition}")
        
        # Create adapted configs for each primer type
        self.adapted_configs = {}
        for primer_type in LAMPPrimerType:
            self.adapted_configs[primer_type] = self.get_adapted_config_for_primer_type(
                primer_type, self.target_composition
            )
        
        return self.adapted_configs
    
    def print_primer_type_summary(self):
        """Print summary of primer-type-specific parameters"""
        
        if not self.adapted_configs:
            print("❌ No adapted configs available")
            return
        
        print(f"\n⚙️ PRIMER-TYPE-SPECIFIC PARAMETERS ({self.target_composition}):")
        print(f"{'='*70}")
        
        # Table header
        print(f"{'Type':<6} {'Tm Range':<12} {'Length':<10} {'GC Range':<12} {'Priority':<8} {'Description'}")
        print(f"{'-'*6} {'-'*12} {'-'*10} {'-'*12} {'-'*8} {'-'*20}")
        
        # Print each primer type
        for primer_type, config in self.adapted_configs.items():
            tm_range = f"{config['tm_range'][0]:.1f}-{config['tm_range'][1]:.1f}°C"
            length_range = f"{config['length_range'][0]}-{config['length_range'][1]}bp"
            gc_range = f"{config['gc_range'][0]:.0f}-{config['gc_range'][1]:.0f}%"
            priority = f"{config['priority_weight']:.1f}"
            description = config['description'][:20] + "..." if len(config['description']) > 20 else config['description']
            
            print(f"{primer_type.value:<6} {tm_range:<12} {length_range:<10} {gc_range:<12} {priority:<8} {description}")
        
        # Print key insights
        print(f"\n💡 KEY INSIGHTS:")
        outer_tm = self.adapted_configs[LAMPPrimerType.F3]['tm_optimal']
        inner_tm = self.adapted_configs[LAMPPrimerType.F1C]['tm_optimal']
        composite_tm = self.adapted_configs[LAMPPrimerType.FIP]['tm_optimal']
        
        print(f"   🌡️  Tm hierarchy: Outer ({outer_tm:.1f}°C) > Composite ({composite_tm:.1f}°C) > Inner ({inner_tm:.1f}°C)")
        print(f"   📏 Length hierarchy: Composite > Outer > Inner")
        print(f"   🎯 Priority: Composite (1.5) > Outer (1.2) > Inner (1.0) > Loop (0.8)")

class DNAThermodynamics:
    """DNA thermodynamics calculations for primer design"""
    
    # Nearest neighbor parameters for enthalpy (kcal/mol)
    NN_ENTHALPY = {
        'AA': -7.9, 'AT': -7.2, 'AC': -8.4, 'AG': -7.8,
        'TA': -7.2, 'TT': -7.9, 'TC': -8.2, 'TG': -8.5,
        'CA': -8.5, 'CT': -7.8, 'CC': -8.0, 'CG': -10.6,
        'GA': -8.2, 'GT': -8.4, 'GC': -9.8, 'GG': -8.0
    }
    
    # Nearest neighbor parameters for entropy (cal/(mol·K))
    NN_ENTROPY = {
        'AA': -22.2, 'AT': -20.4, 'AC': -22.4, 'AG': -21.0,
        'TA': -21.3, 'TT': -22.2, 'TC': -22.2, 'TG': -22.7,
        'CA': -22.7, 'CT': -21.0, 'CC': -19.9, 'CG': -27.2,
        'GA': -22.2, 'GT': -22.4, 'GC': -24.4, 'GG': -19.9
    }
    
    # Nearest neighbor parameters for Gibbs free energy at 37°C (kcal/mol)
    NN_GIBBS_37 = {
        'AA': -1.00, 'AT': -0.88, 'AC': -1.45, 'AG': -1.28,
        'TA': -0.58, 'TT': -1.00, 'TC': -1.30, 'TG': -1.45,
        'CA': -1.45, 'CT': -1.28, 'CC': -1.84, 'CG': -2.17,
        'GA': -1.30, 'GT': -1.45, 'GC': -2.24, 'GG': -1.84
    }
    
    @staticmethod
    def calculate_tm(sequence: str, na_conc: float = 50.0, mg_conc: float = 8.0) -> float:
        """Calculate melting temperature using nearest neighbor method"""
        if len(sequence) < 2:
            return 0.0
        
        sequence = sequence.upper()
        
        # Calculate enthalpy and entropy
        dh = 0.0  # kcal/mol
        ds = 0.0  # cal/(mol·K)
        
        # Add initiation parameters
        if sequence[0] in 'AT':
            dh += 2.3
            ds += 4.1
        else:  # GC
            dh += 0.1
            ds -= 2.8
        
        if sequence[-1] in 'AT':
            dh += 2.3
            ds += 4.1
        else:  # GC
            dh += 0.1
            ds -= 2.8
        
        # Add nearest neighbor contributions
        for i in range(len(sequence) - 1):
            nn = sequence[i:i+2]
            if nn in DNAThermodynamics.NN_ENTHALPY:
                dh += DNAThermodynamics.NN_ENTHALPY[nn]
                ds += DNAThermodynamics.NN_ENTROPY[nn]
        
        # Salt correction (simplified)
        salt_effect = 0.368 * math.log(na_conc / 1000.0)
        
        # Calculate Tm
        tm = (dh * 1000) / (ds + 1.987 * math.log(0.25e-6)) - 273.15 + salt_effect
        
        return tm
    
    @staticmethod
    def calculate_gibbs_free_energy(sequence: str, temperature: float = 37.0) -> float:
        """Calculate Gibbs free energy of DNA duplex formation (kcal/mol)"""
        if len(sequence) < 2:
            return 0.0
        
        sequence = sequence.upper()
        dg = 0.0
        
        # Initiation penalty
        dg += 0.98  # kcal/mol for initiation
        
        # Add nearest neighbor contributions
        for i in range(len(sequence) - 1):
            nn = sequence[i:i+2]
            if nn in DNAThermodynamics.NN_GIBBS_37:
                dg += DNAThermodynamics.NN_GIBBS_37[nn]
        
        # Temperature correction (simplified)
        if temperature != 37.0:
            # Basic temperature correction using Van't Hoff equation
            temp_correction = (temperature - 37.0) * 0.01  # Approximate
            dg += temp_correction
        
        return dg
    
    @staticmethod
    def calculate_3prime_stability(sequence: str, length: int = 5) -> float:
        """Calculate stability of 3' end (last 'length' bases)"""
        if len(sequence) < length:
            return DNAThermodynamics.calculate_gibbs_free_energy(sequence)
        
        three_prime = sequence[-length:]
        return DNAThermodynamics.calculate_gibbs_free_energy(three_prime)
    
    @staticmethod
    def calculate_hairpin_dg(sequence: str) -> float:
        """Calculate minimum Gibbs free energy for hairpin formation"""
        min_dg = 0.0
        sequence = sequence.upper()
        
        # Check for potential hairpin structures
        # Simplified approach: look for inverted repeats
        for i in range(len(sequence) - 6):  # Minimum hairpin size
            for j in range(i + 6, len(sequence)):
                # Check for complementary regions
                stem1 = sequence[i:i+3]
                stem2 = sequence[j-2:j+1]
                
                if DNAThermodynamics._is_complementary(stem1, stem2):
                    # Calculate stability of this potential hairpin
                    stem_dg = DNAThermodynamics.calculate_gibbs_free_energy(stem1)
                    loop_penalty = 3.0  # kcal/mol penalty for loop formation
                    hairpin_dg = stem_dg - loop_penalty
                    
                    if hairpin_dg < min_dg:
                        min_dg = hairpin_dg
        
        return min_dg
    
    @staticmethod
    def calculate_dimer_dg(seq1: str, seq2: str) -> float:
        """Calculate Gibbs free energy for primer-dimer formation"""
        min_dg = 0.0
        
        # Check all possible alignments between primers
        for offset in range(-len(seq2) + 3, len(seq1) - 2):
            dg = DNAThermodynamics._calculate_alignment_dg(seq1, seq2, offset)
            if dg < min_dg:
                min_dg = dg
        
        return min_dg
    
    @staticmethod
    def _is_complementary(seq1: str, seq2: str) -> bool:
        """Check if two sequences are complementary"""
        complement = {'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G'}
        if len(seq1) != len(seq2):
            return False
        
        for i in range(len(seq1)):
            if complement.get(seq1[i]) != seq2[len(seq2)-1-i]:
                return False
        return True
    
    @staticmethod
    def _calculate_alignment_dg(seq1: str, seq2: str, offset: int) -> float:
        """Calculate dG for a specific alignment between two sequences"""
        dg = 0.0
        matches = 0
        
        for i in range(len(seq1)):
            j = i - offset
            if 0 <= j < len(seq2):
                if DNAThermodynamics._bases_complement(seq1[i], seq2[j]):
                    matches += 1
                    # Add stabilizing energy for each match
                    dg -= 1.0  # Simplified: -1 kcal/mol per match
        
        # Penalty for mismatches and gaps
        if matches < 3:  # Minimum stable interaction
            return 0.0
        
        # Additional penalty for end effects
        dg += 2.0  # Initiation penalty
        
        return dg
    
    @staticmethod
    def _bases_complement(base1: str, base2: str) -> bool:
        """Check if two bases are complementary"""
        complements = {'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G'}
        return complements.get(base1.upper()) == base2.upper()
    
    @staticmethod
    def calculate_gc_content(sequence: str) -> float:
        """Calculate GC content percentage"""
        sequence = sequence.upper()
        gc_count = sequence.count('G') + sequence.count('C')
        return (gc_count / len(sequence)) * 100.0 if sequence else 0.0

class SequenceProcessor:
    """Process sequences for different sequencing approaches"""
    
    @staticmethod
    def process_bisulfite(sequence: str) -> str:
        """Convert sequence for bisulfite sequencing (C -> T, except CpG)"""
        sequence = sequence.upper()
        processed = list(sequence)
        
        for i, base in enumerate(sequence):
            if base == 'C':
                # Check if it's part of CpG
                if i < len(sequence) - 1 and sequence[i + 1] == 'G':
                    continue  # Keep CpG methylation sites
                else:
                    processed[i] = 'T'  # Convert non-CpG cytosines
        
        return ''.join(processed)
    
    @staticmethod
    def reverse_complement(sequence: str) -> str:
        """Return reverse complement of DNA sequence"""
        complement = {'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G',
                     'N': 'N', 'R': 'Y', 'Y': 'R', 'S': 'S',
                     'W': 'W', 'K': 'M', 'M': 'K'}
        
        return ''.join(complement.get(base.upper(), base) 
                      for base in reversed(sequence))

class LoopPrimerDesigner:
    """Enhanced designer for LAMP loop primers with proper topology understanding"""
    
    def __init__(self, config: PrimerConfig):
        self.config = config
        self.thermo = DNAThermodynamics()
        self.processor = SequenceProcessor()
        
    def design_loop_primers(self, target_sequence: str, lamp_set: LAMPPrimerSet, 
                          debug: bool = False) -> Tuple[Optional[Primer], Optional[Primer]]:
        """Design loop primers with proper understanding of LAMP structure"""
        
        if debug:
            print(f"\n🔄 ENHANCED LOOP PRIMER DESIGN")
            print(f"{'='*50}")
        
        # Extract component positions from composite primers
        # This is approximate since we don't store original component positions
        f3_pos = lamp_set.f3.start
        b3_pos = lamp_set.b3.start
        
        # Estimate F2 and F1 positions from FIP
        # FIP = F1c + F2, so F2 is at the 3' end of FIP region
        fip_start = lamp_set.fip.start
        fip_end = lamp_set.fip.end
        
        # Estimate B2 and B1 positions from BIP
        bip_start = lamp_set.bip.start
        bip_end = lamp_set.bip.end
        
        if debug:
            print(f"📍 Primer positions:")
            print(f"   F3: {f3_pos}-{lamp_set.f3.end}")
            print(f"   FIP region: {fip_start}-{fip_end}")
            print(f"   BIP region: {bip_start}-{bip_end}")
            print(f"   B3: {b3_pos}-{lamp_set.b3.end}")
        
        # Design LF (binds between F2 and F1)
        lf_primer = self._design_lf_primer(target_sequence, lamp_set, debug)
        
        # Design LB (binds between B1 and B2)
        lb_primer = self._design_lb_primer(target_sequence, lamp_set, debug)
        
        if debug:
            if lf_primer:
                print(f"\n✅ LF designed: {lf_primer.sequence} (Tm={lf_primer.tm:.1f}°C)")
            else:
                print(f"\n❌ LF design failed")
                
            if lb_primer:
                print(f"✅ LB designed: {lb_primer.sequence} (Tm={lb_primer.tm:.1f}°C)")
            else:
                print(f"❌ LB design failed")
        
        return lf_primer, lb_primer
    
    def _design_lf_primer(self, sequence: str, lamp_set: LAMPPrimerSet, debug: bool) -> Optional[Primer]:
        """Design LF primer for the forward loop region"""
        
        # LF binds to the loop between F2 and F1
        # This is typically 20-60bp downstream of F2
        
        # Estimate F2 end position (F2 is part of FIP, at the 3' end)
        # Since FIP = F1c + F2, and F2 binds near the F3 region
        # We look for LF binding sites between F3 end and FIP start
        
        search_start = lamp_set.f3.end + 10  # At least 10bp gap from F3
        search_end = lamp_set.fip.start - 10  # Before FIP binding region
        
        if search_end <= search_start:
            # Try alternative: look downstream of FIP
            search_start = lamp_set.fip.end + 5
            search_end = min(search_start + 100, lamp_set.b3.start - 20)
        
        if debug:
            print(f"\n🔍 LF search region: {search_start}-{search_end}")
            print(f"   Region length: {search_end - search_start}bp")
        
        if search_end - search_start < self.config.loop_primer_length:
            if debug:
                print(f"   ❌ Search region too small for loop primers")
            return None
        
        # Find best loop primer in this region
        candidates = self._find_loop_candidates_in_region(
            sequence, search_start, search_end, "LF", debug
        )
        
        if not candidates:
            # Try with relaxed constraints
            candidates = self._find_loop_candidates_relaxed(
                sequence, search_start, search_end, "LF", debug
            )
        
        if candidates:
            # Select best based on Tm compatibility with main primers
            best_candidate = self._select_best_loop_candidate(candidates, lamp_set, "LF")
            return best_candidate
        
        return None
    
    def _design_lb_primer(self, sequence: str, lamp_set: LAMPPrimerSet, debug: bool) -> Optional[Primer]:
        """Design LB primer for the backward loop region"""
        
        # LB binds to the loop between B1 and B2
        # This is typically in the region between BIP and B3
        
        search_start = lamp_set.bip.end + 10
        search_end = lamp_set.b3.start - 10
        
        if search_end <= search_start:
            # Try alternative: look upstream of BIP
            search_end = lamp_set.bip.start - 5
            search_start = max(lamp_set.fip.end + 20, search_end - 100)
        
        if debug:
            print(f"\n🔍 LB search region: {search_start}-{search_end}")
            print(f"   Region length: {search_end - search_start}bp")
        
        if search_end - search_start < self.config.loop_primer_length:
            if debug:
                print(f"   ❌ Search region too small for loop primers")
            return None
        
        # Find best loop primer in this region
        candidates = self._find_loop_candidates_in_region(
            sequence, search_start, search_end, "LB", debug
        )
        
        if not candidates:
            # Try with relaxed constraints
            candidates = self._find_loop_candidates_relaxed(
                sequence, search_start, search_end, "LB", debug
            )
        
        if candidates:
            # Select best based on Tm compatibility
            best_candidate = self._select_best_loop_candidate(candidates, lamp_set, "LB")
            return best_candidate
        
        return None
    
    def _find_loop_candidates_in_region(self, sequence: str, start_pos: int, 
                                      end_pos: int, primer_type: str, 
                                      debug: bool) -> List[Primer]:
        """Find loop primer candidates in a specific region"""
        
        candidates = []
        
        # Try different primer lengths
        for length in range(self.config.loop_primer_length - 3, 
                          self.config.loop_primer_length + 4):
            
            # Scan the region
            for pos in range(start_pos, min(end_pos - length + 1, len(sequence) - length + 1)):
                primer_seq = sequence[pos:pos + length]
                
                # Skip if contains N
                if 'N' in primer_seq:
                    continue
                
                # Calculate properties
                tm = self.thermo.calculate_tm(primer_seq, self.config.na_conc, self.config.mg_conc)
                gc = self.thermo.calculate_gc_content(primer_seq)
                dg_3prime = self.thermo.calculate_3prime_stability(primer_seq)
                dg_hairpin = self.thermo.calculate_hairpin_dg(primer_seq)
                
                # Loop primers can have more relaxed constraints
                tm_ok = self.config.min_tm - 10 <= tm <= self.config.max_tm + 5
                gc_ok = 20 <= gc <= 80  # Very relaxed GC
                quality_ok = self._check_loop_primer_quality(primer_seq)
                
                if tm_ok and gc_ok and quality_ok:
                    primer = Primer(
                        sequence=primer_seq,
                        start=pos,
                        end=pos + length,
                        tm=tm,
                        gc_content=gc,
                        primer_type=primer_type,
                        dg_3prime=dg_3prime,
                        dg_hairpin=dg_hairpin
                    )
                    candidates.append(primer)
        
        if debug:
            print(f"   Found {len(candidates)} candidates in region")
            if candidates:
                tms = [c.tm for c in candidates]
                print(f"   Tm range: {min(tms):.1f}-{max(tms):.1f}°C")
        
        return candidates
    
    def _find_loop_candidates_relaxed(self, sequence: str, start_pos: int,
                                    end_pos: int, primer_type: str,
                                    debug: bool) -> List[Primer]:
        """Find loop candidates with very relaxed constraints"""
        
        if debug:
            print(f"   🔄 Trying relaxed constraints...")
        
        candidates = []
        
        # Very relaxed length range
        for length in range(15, 26):
            for pos in range(start_pos, min(end_pos - length + 1, len(sequence) - length + 1)):
                primer_seq = sequence[pos:pos + length]
                
                if 'N' in primer_seq:
                    continue
                
                # Super relaxed - just check basic viability
                tm = self.thermo.calculate_tm(primer_seq)
                gc = self.thermo.calculate_gc_content(primer_seq)
                
                # Very relaxed constraints
                if 40 <= tm <= 85 and 15 <= gc <= 85:
                    # Avoid really bad sequences
                    if not re.search(r'(.)\1{5,}', primer_seq):  # No 6+ homopolymers
                        primer = Primer(
                            sequence=primer_seq,
                            start=pos,
                            end=pos + length,
                            tm=tm,
                            gc_content=gc,
                            primer_type=primer_type,
                            dg_3prime=0,
                            dg_hairpin=0
                        )
                        candidates.append(primer)
        
        if debug and candidates:
            print(f"   Found {len(candidates)} candidates with relaxed constraints")
        
        return candidates
    
    def _check_loop_primer_quality(self, sequence: str) -> bool:
        """Check loop primer quality with appropriate constraints"""
        
        # More relaxed than main primers
        if re.search(r'(.)\1{5,}', sequence):  # No 6+ homopolymers
            return False
        
        # Check for extreme GC skew in small windows
        if len(sequence) >= 10:
            for i in range(len(sequence) - 9):
                window = sequence[i:i+10]
                window_gc = (window.count('G') + window.count('C')) / 10
                if window_gc == 0 or window_gc == 1:  # All AT or all GC
                    return False
        
        return True
    
    def _select_best_loop_candidate(self, candidates: List[Primer], 
                                  lamp_set: LAMPPrimerSet,
                                  loop_type: str) -> Optional[Primer]:
        """Select best loop primer based on compatibility"""
        
        if not candidates:
            return None
        
        # Get main primer Tms for comparison
        main_tms = [lamp_set.f3.tm, lamp_set.b3.tm, lamp_set.fip.tm, lamp_set.bip.tm]
        avg_tm = sum(main_tms) / len(main_tms)
        
        # Score candidates
        scored = []
        for candidate in candidates:
            # Tm compatibility (prefer similar to main primers)
            tm_diff = abs(candidate.tm - avg_tm)
            tm_score = 1.0 / (1.0 + tm_diff / 5.0)
            
            # Position score (prefer middle of search region)
            if loop_type == "LF":
                # For LF, prefer candidates not too close to F3 or FIP
                pos_score = 1.0
                if abs(candidate.start - lamp_set.f3.end) < 15:
                    pos_score *= 0.7
                if abs(candidate.end - lamp_set.fip.start) < 15:
                    pos_score *= 0.7
            else:  # LB
                # For LB, prefer candidates not too close to BIP or B3
                pos_score = 1.0
                if abs(candidate.start - lamp_set.bip.end) < 15:
                    pos_score *= 0.7
                if abs(candidate.end - lamp_set.b3.start) < 15:
                    pos_score *= 0.7
            
            # GC balance score
            gc_score = 1.0
            if candidate.gc_content < 30 or candidate.gc_content > 70:
                gc_score = 0.8
            
            # Overall score
            total_score = (0.5 * tm_score + 0.3 * pos_score + 0.2 * gc_score)
            scored.append((candidate, total_score))
        
        # Sort by score and return best
        scored.sort(key=lambda x: x[1], reverse=True)
        return scored[0][0]
    
    def analyze_loop_primer_failure(self, target_sequence: str, lamp_set: LAMPPrimerSet):
        """Analyze why loop primer design failed"""
        
        print(f"\n🔍 LOOP PRIMER FAILURE ANALYSIS")
        print(f"{'='*50}")
        
        # Analyze primer positions
        print(f"\n📍 Main primer positions:")
        print(f"   F3: {lamp_set.f3.start}-{lamp_set.f3.end} ({lamp_set.f3.end - lamp_set.f3.start}bp)")
        print(f"   FIP: {lamp_set.fip.start}-{lamp_set.fip.end} ({lamp_set.fip.end - lamp_set.fip.start}bp)")
        print(f"   BIP: {lamp_set.bip.start}-{lamp_set.bip.end} ({lamp_set.bip.end - lamp_set.bip.start}bp)")
        print(f"   B3: {lamp_set.b3.start}-{lamp_set.b3.end} ({lamp_set.b3.end - lamp_set.b3.start}bp)")
        
        # Calculate gaps
        f3_fip_gap = lamp_set.fip.start - lamp_set.f3.end
        fip_bip_gap = lamp_set.bip.start - lamp_set.fip.end
        bip_b3_gap = lamp_set.b3.start - lamp_set.bip.end
        
        print(f"\n📏 Inter-primer gaps:")
        print(f"   F3-FIP gap: {f3_fip_gap}bp")
        print(f"   FIP-BIP gap: {fip_bip_gap}bp")
        print(f"   BIP-B3 gap: {bip_b3_gap}bp")
        
        # Analyze potential loop regions
        print(f"\n🔄 Potential loop regions:")
        
        # LF region analysis
        if f3_fip_gap > 20:
            lf_region = target_sequence[lamp_set.f3.end:lamp_set.fip.start]
            print(f"\n   LF region (F3-FIP gap): {len(lf_region)}bp")
            print(f"   GC content: {self.thermo.calculate_gc_content(lf_region):.1f}%")
            print(f"   Sequence preview: {lf_region[:30]}...")
            
            # Test sample primers
            if len(lf_region) >= 20:
                test_seq = lf_region[:20]
                test_tm = self.thermo.calculate_tm(test_seq)
                print(f"   Sample 20bp Tm: {test_tm:.1f}°C")
        else:
            print(f"   ❌ F3-FIP gap too small for LF ({f3_fip_gap}bp)")
        
        # LB region analysis  
        if bip_b3_gap > 20:
            lb_region = target_sequence[lamp_set.bip.end:lamp_set.b3.start]
            print(f"\n   LB region (BIP-B3 gap): {len(lb_region)}bp")
            print(f"   GC content: {self.thermo.calculate_gc_content(lb_region):.1f}%")
            print(f"   Sequence preview: {lb_region[:30]}...")
            
            # Test sample primers
            if len(lb_region) >= 20:
                test_seq = lb_region[:20]
                test_tm = self.thermo.calculate_tm(test_seq)
                print(f"   Sample 20bp Tm: {test_tm:.1f}°C")
        else:
            print(f"   ❌ BIP-B3 gap too small for LB ({bip_b3_gap}bp)")
        
        print(f"\n💡 Suggestions:")
        print(f"   - Ensure main primers have sufficient gaps (>30bp)")
        print(f"   - Consider adjusting main primer positions")
        print(f"   - Try manual loop primer design in identified regions")

@dataclass
class OffTargetMatch:
    """Represents an off-target match for a primer"""
    sequence_id: str
    position: int
    primer_sequence: str
    target_sequence: str
    tm: float
    mismatches: int
    gaps: int
    alignment_score: float
    
@dataclass
class SpecificityResult:
    """Results of specificity analysis"""
    primer: Primer
    off_targets: List[OffTargetMatch]
    is_specific: bool
    max_off_target_tm: float
    total_off_targets: int

@dataclass
class DesignResults:
    """Container for primer design results - FIXED VERSION with optional fields"""
    primer_sets: List[LAMPPrimerSet]
    total_candidates: int
    total_combinations: int
    filtering_stats: Dict[str, int]
    scores: List[float]
    specificity_results: Optional[Dict[str, List[SpecificityResult]]] = None
    
    # NEW: Optional type-specific attributes (fixes the constructor error)
    composition_analysis: Optional[Dict] = None
    design_method: Optional[str] = None
    
    def get_top_n(self, n: int = 3) -> List[LAMPPrimerSet]:
        """Get top N primer sets"""
        return self.primer_sets[:n]
    
    def get_all(self) -> List[LAMPPrimerSet]:
        """Get all valid primer sets"""
        return self.primer_sets
    
    def get_statistics(self) -> Dict:
        """Get design statistics"""
        stats = {
            'total_candidates': self.total_candidates,
            'total_combinations': self.total_combinations,
            'valid_sets': len(self.primer_sets),
            'filtering_stats': self.filtering_stats,
            'score_range': (min(self.scores) if self.scores else 0, 
                          max(self.scores) if self.scores else 0)
        }
        
        if self.specificity_results:
            stats['specificity_stats'] = {
                primer_type: {
                    'total_checked': len(results),
                    'specific_primers': sum(1 for r in results if r.is_specific),
                    'avg_off_targets': np.mean([r.total_off_targets for r in results]) if results else 0
                }
                for primer_type, results in self.specificity_results.items()
            }
        
        # Add type-specific stats if available
        if self.design_method:
            stats['design_method'] = self.design_method
        if self.composition_analysis:
            stats['composition_analysis'] = self.composition_analysis
        
        return stats
    
@dataclass
class BatchTarget:
    """Represents a target sequence for batch processing"""
    target_id: str
    sequence: str
    description: Optional[str] = None
    seq_type: Optional[SequenceType] = None
    custom_config: Optional[Dict] = None

@dataclass
class BatchResults:
    """Results from batch primer design"""
    target_results: Dict[str, DesignResults]
    batch_statistics: Dict[str, any]
    cross_target_analysis: Dict[str, any]
    batch_effects: Dict[str, any]
    
    def get_best_primers_per_target(self) -> Dict[str, LAMPPrimerSet]:
        """Get best primer set for each target"""
        best_primers = {}
        for target_id, results in self.target_results.items():
            if results.primer_sets:
                best_primers[target_id] = results.primer_sets[0]
        return best_primers
    
    def get_batch_summary(self) -> Dict:
        """Get summary statistics across all targets"""
        total_sets = sum(len(r.primer_sets) for r in self.target_results.values())
        successful_targets = sum(1 for r in self.target_results.values() if r.primer_sets)
        
        return {
            'total_targets': len(self.target_results),
            'successful_targets': successful_targets,
            'total_primer_sets': total_sets,
            'success_rate': successful_targets / len(self.target_results) if self.target_results else 0,
            'batch_effects': self.batch_effects,
            'cross_target_interactions': self.cross_target_analysis
        }

class BatchPrimerDesigner:
    """Handles batch processing of multiple targets with batch effect analysis"""
    
    def __init__(self, config: PrimerConfig):
        self.config = config
        self.individual_designer = LAMPPrimerDesigner(config)
        self.loop_designer = LoopPrimerDesigner(config)
        
    def design_batch_primers(self, targets: List[BatchTarget], 
                           include_loop_primers: bool = True,
                           reference_sequences: Optional[Dict[str, str]] = None) -> BatchResults:
        """Design primers for multiple targets with batch effect analysis"""
        
        print(f"🧬 Starting batch primer design for {len(targets)} targets...")
        
        target_results = {}
        all_primer_sets = []
        
        # Design primers for each target
        for i, target in enumerate(targets):
            print(f"Processing target {i+1}/{len(targets)}: {target.target_id}")
            
            # Use target-specific config if provided
            if target.custom_config:
                target_config = self._merge_configs(self.config, target.custom_config)
                target_designer = LAMPPrimerDesigner(target_config)
            else:
                target_designer = self.individual_designer
            
            # Override sequence type if specified
            if target.seq_type:
                target_designer.config.seq_type = target.seq_type
            
            # Design primers
            results = target_designer.design_primers(
                target.sequence,
                ResultFormat.DETAILED,
                reference_sequences=reference_sequences
            )
            
            # Add loop primers if requested and primers were found
            if include_loop_primers and isinstance(results, DesignResults) and results.primer_sets:
                enhanced_sets = []
                for lamp_set in results.primer_sets[:5]:  # Add loops to top 5 sets
                    lf, lb = self.loop_designer.design_loop_primers(target.sequence, lamp_set)
                    enhanced_set = LAMPPrimerSet(
                        lamp_set.f3, lamp_set.b3, lamp_set.fip, lamp_set.bip, lf, lb
                    )
                    enhanced_sets.append(enhanced_set)
                
                # Replace original sets with enhanced ones
                if enhanced_sets:
                    results.primer_sets[:len(enhanced_sets)] = enhanced_sets
            
            target_results[target.target_id] = results
            
            # Collect all primer sets for batch analysis
            if isinstance(results, DesignResults):
                all_primer_sets.extend(results.primer_sets)
        
        # Perform batch effect analysis
        print("📊 Analyzing batch effects...")
        batch_effects = self._analyze_batch_effects(target_results)
        
        # Cross-target interaction analysis
        print("🔍 Checking cross-target interactions...")
        cross_target_analysis = self._analyze_cross_target_interactions(target_results)
        
        # Batch statistics
        batch_stats = self._calculate_batch_statistics(target_results)
        
        return BatchResults(
            target_results=target_results,
            batch_statistics=batch_stats,
            cross_target_analysis=cross_target_analysis,
            batch_effects=batch_effects
        )
    
    def _merge_configs(self, base_config: PrimerConfig, custom_config: Dict) -> PrimerConfig:
        """Merge custom configuration with base configuration"""
        config_dict = base_config.__dict__.copy()
        config_dict.update(custom_config)
        return PrimerConfig(**config_dict)
    
    def _analyze_batch_effects(self, target_results: Dict[str, DesignResults]) -> Dict[str, any]:
        """Analyze potential batch effects across targets"""
        
        all_scores = []
        all_tms = []
        all_gcs = []
        target_success_rates = []
        
        for target_id, results in target_results.items():
            if isinstance(results, DesignResults) and results.primer_sets:
                # Collect metrics
                target_scores = results.scores if results.scores else []
                all_scores.extend(target_scores)
                
                # Collect Tm and GC values
                for primer_set in results.primer_sets:
                    primers = [primer_set.f3, primer_set.b3, primer_set.fip, primer_set.bip]
                    all_tms.extend([p.tm for p in primers])
                    all_gcs.extend([p.gc_content for p in primers])
                
                target_success_rates.append(len(results.primer_sets))
            else:
                target_success_rates.append(0)
        
        # Calculate batch effect metrics
        batch_effects = {
            'score_variation': {
                'mean': np.mean(all_scores) if all_scores else 0,
                'std': np.std(all_scores) if all_scores else 0,
                'cv': np.std(all_scores) / np.mean(all_scores) if all_scores and np.mean(all_scores) > 0 else 0
            },
            'tm_variation': {
                'mean': np.mean(all_tms) if all_tms else 0,
                'std': np.std(all_tms) if all_tms else 0,
                'range': (min(all_tms), max(all_tms)) if all_tms else (0, 0)
            },
            'gc_variation': {
                'mean': np.mean(all_gcs) if all_gcs else 0,
                'std': np.std(all_gcs) if all_gcs else 0,
                'range': (min(all_gcs), max(all_gcs)) if all_gcs else (0, 0)
            },
            'success_rate_variation': {
                'mean': np.mean(target_success_rates),
                'std': np.std(target_success_rates),
                'min_max': (min(target_success_rates), max(target_success_rates))
            }
        }
        
        return batch_effects
    
    def _analyze_cross_target_interactions(self, target_results: Dict[str, DesignResults]) -> Dict[str, any]:
        """Analyze potential cross-reactions between primers from different targets"""
        
        interactions = []
        target_primers = {}
        
        # Collect best primer set from each target
        for target_id, results in target_results.items():
            if isinstance(results, DesignResults) and results.primer_sets:
                best_set = results.primer_sets[0]
                target_primers[target_id] = [best_set.f3, best_set.b3, best_set.fip, best_set.bip]
        
        # Check cross-interactions
        target_ids = list(target_primers.keys())
        for i, target1 in enumerate(target_ids):
            for j, target2 in enumerate(target_ids[i+1:], i+1):
                
                primers1 = target_primers[target1]
                primers2 = target_primers[target2]
                
                # Check all primer combinations between targets
                for primer1 in primers1:
                    for primer2 in primers2:
                        dimer_dg = DNAThermodynamics.calculate_dimer_dg(
                            primer1.sequence, primer2.sequence
                        )
                        
                        if dimer_dg < self.config.max_dimer_dg:
                            interactions.append({
                                'target1': target1,
                                'target2': target2,
                                'primer1': f"{primer1.primer_type}:{primer1.sequence}",
                                'primer2': f"{primer2.primer_type}:{primer2.sequence}",
                                'dimer_dg': dimer_dg
                            })
        
        return {
            'total_interactions': len(interactions),
            'problematic_pairs': len([i for i in interactions if i['dimer_dg'] < -7.0]),
            'interactions': interactions[:10]  # Keep top 10 for review
        }
    
    def _calculate_batch_statistics(self, target_results: Dict[str, DesignResults]) -> Dict[str, any]:
        """Calculate comprehensive batch statistics"""
        
        total_targets = len(target_results)
        successful_targets = 0
        total_primer_sets = 0
        total_candidates = 0
        
        filtering_summary = defaultdict(int)
        
        for results in target_results.values():
            if isinstance(results, DesignResults):
                if results.primer_sets:
                    successful_targets += 1
                    total_primer_sets += len(results.primer_sets)
                
                total_candidates += results.total_candidates
                
                # Aggregate filtering statistics
                for filter_type, count in results.filtering_stats.items():
                    filtering_summary[filter_type] += count
        
        return {
            'total_targets': total_targets,
            'successful_targets': successful_targets,
            'success_rate': successful_targets / total_targets if total_targets > 0 else 0,
            'total_primer_sets': total_primer_sets,
            'avg_sets_per_target': total_primer_sets / successful_targets if successful_targets > 0 else 0,
            'total_candidates_evaluated': total_candidates,
            'filtering_summary': dict(filtering_summary)
        }

class NonSpecificityChecker:
    """Checks primers for non-specific binding using k-mer indexing and alignment"""
    
    def __init__(self, config: PrimerConfig):
        self.config = config
        self.kmer_index: Dict[str, List[Tuple[str, int]]] = defaultdict(list)
        self.reference_sequences: Dict[str, str] = {}
        self.thermo = DNAThermodynamics()
    
    def build_kmer_index(self, reference_sequences: Dict[str, str]):
        """Build k-mer index from reference sequences for fast screening"""
        self.reference_sequences = reference_sequences
        self.kmer_index.clear()
        
        for seq_id, sequence in reference_sequences.items():
            sequence = sequence.upper()
            for i in range(len(sequence) - self.config.kmer_size + 1):
                kmer = sequence[i:i + self.config.kmer_size]
                self.kmer_index[kmer].append((seq_id, i))
    
    def find_potential_off_targets(self, primer_sequence: str) -> List[Tuple[str, int]]:
        """Find potential off-target sites using k-mer screening"""
        primer_seq = primer_sequence.upper()
        potential_sites = set()
        
        # Check all k-mers in the primer
        for i in range(len(primer_seq) - self.config.kmer_size + 1):
            kmer = primer_seq[i:i + self.config.kmer_size]
            if kmer in self.kmer_index:
                potential_sites.update(self.kmer_index[kmer])
        
        return list(potential_sites)
    
    def align_primer_to_target(self, primer_seq: str, target_seq: str, 
                              target_pos: int) -> Optional[OffTargetMatch]:
        """Perform local alignment between primer and potential target site"""
        primer_len = len(primer_seq)
        
        # Check region around the k-mer match
        start_pos = max(0, target_pos - 5)
        end_pos = min(len(target_seq), target_pos + primer_len + 5)
        target_region = target_seq[start_pos:end_pos]
        
        best_match = None
        best_score = -float('inf')
        
        # Try different alignments within the region
        for offset in range(len(target_region) - primer_len + 1):
            target_subseq = target_region[offset:offset + primer_len]
            if len(target_subseq) != primer_len:
                continue
            
            # Calculate alignment score
            matches, mismatches, gaps = self._calculate_alignment_stats(primer_seq, target_subseq)
            
            # Simple scoring: +2 for match, -1 for mismatch
            score = matches * 2 - mismatches * 1
            
            if score > best_score and matches >= self.config.kmer_size:
                # Calculate Tm for this alignment
                tm = self._calculate_alignment_tm(primer_seq, target_subseq)
                
                if tm >= self.config.min_off_target_tm:
                    best_match = OffTargetMatch(
                        sequence_id="",  # Will be filled by caller
                        position=start_pos + offset,
                        primer_sequence=primer_seq,
                        target_sequence=target_subseq,
                        tm=tm,
                        mismatches=mismatches,
                        gaps=gaps,
                        alignment_score=score
                    )
                    best_score = score
        
        return best_match
    
    def _calculate_alignment_stats(self, seq1: str, seq2: str) -> Tuple[int, int, int]:
        """Calculate matches, mismatches, and gaps between two sequences"""
        matches = mismatches = gaps = 0
        
        for i in range(min(len(seq1), len(seq2))):
            if seq1[i] == seq2[i]:
                matches += 1
            else:
                mismatches += 1
        
        gaps = abs(len(seq1) - len(seq2))
        return matches, mismatches, gaps
    
    def _calculate_alignment_tm(self, primer_seq: str, target_seq: str) -> float:
        """Calculate melting temperature for primer-target alignment"""
        # For mismatched sequences, reduce Tm based on mismatches
        matches, mismatches, gaps = self._calculate_alignment_stats(primer_seq, target_seq)
        
        if matches == 0:
            return 0.0
        
        # Calculate base Tm for matched regions
        matched_portion = matches / len(primer_seq)
        base_tm = self.thermo.calculate_tm(primer_seq)
        
        # Penalty for mismatches and gaps
        mismatch_penalty = mismatches * 5.0  # 5°C penalty per mismatch
        gap_penalty = gaps * 3.0  # 3°C penalty per gap
        
        adjusted_tm = base_tm * matched_portion - mismatch_penalty - gap_penalty
        
        return max(0.0, adjusted_tm)
    
    def check_primer_specificity(self, primer: Primer) -> SpecificityResult:
        """Check specificity of a single primer"""
        if not self.reference_sequences:
            # No reference database loaded
            return SpecificityResult(
                primer=primer,
                off_targets=[],
                is_specific=True,
                max_off_target_tm=0.0,
                total_off_targets=0
            )
        
        # Find potential off-target sites
        potential_sites = self.find_potential_off_targets(primer.sequence)
        
        off_targets = []
        max_off_target_tm = 0.0
        
        # Analyze each potential site
        for seq_id, pos in potential_sites:
            target_seq = self.reference_sequences[seq_id]
            match = self.align_primer_to_target(primer.sequence, target_seq, pos)
            
            if match:
                match.sequence_id = seq_id
                off_targets.append(match)
                max_off_target_tm = max(max_off_target_tm, match.tm)
        
        # Sort by Tm (highest first)
        off_targets.sort(key=lambda x: x.tm, reverse=True)
        
        # Determine specificity
        significant_off_targets = [ot for ot in off_targets 
                                 if ot.tm >= self.config.min_off_target_tm]
        
        is_specific = len(significant_off_targets) <= self.config.max_off_targets
        
        return SpecificityResult(
            primer=primer,
            off_targets=off_targets[:10],  # Keep top 10 for analysis
            is_specific=is_specific,
            max_off_target_tm=max_off_target_tm,
            total_off_targets=len(significant_off_targets)
        )
    
    def check_primer_set_specificity(self, lamp_set: LAMPPrimerSet) -> Dict[str, SpecificityResult]:
        """Check specificity for entire LAMP primer set"""
        primers_to_check = {
            'F3': lamp_set.f3,
            'B3': lamp_set.b3,
            'FIP': lamp_set.fip,
            'BIP': lamp_set.bip
        }
        
        if lamp_set.lf:
            primers_to_check['LF'] = lamp_set.lf
        if lamp_set.lb:
            primers_to_check['LB'] = lamp_set.lb
        
        results = {}
        for primer_type, primer in primers_to_check.items():
            results[primer_type] = self.check_primer_specificity(primer)
        
        return results
    
    def load_reference_database(self, database_path: Optional[str] = None, 
                               sequences: Optional[Dict[str, str]] = None):
        """Load reference sequences for specificity checking"""
        if sequences:
            # Use provided sequences
            self.build_kmer_index(sequences)
        elif database_path:
            # Load from file (FASTA format)
            try:
                sequences = self._load_fasta_file(database_path)
                self.build_kmer_index(sequences)
                print(f"Loaded {len(sequences)} reference sequences from {database_path}")
            except Exception as e:
                print(f"Error loading database: {e}")
        else:
            print("No reference database provided. Specificity checking disabled.")
    
    def _load_fasta_file(self, filepath: str) -> Dict[str, str]:
        """Load sequences from FASTA file"""
        sequences = {}
        current_id = None
        current_seq = []
        
        try:
            with open(filepath, 'r') as f:
                for line in f:
                    line = line.strip()
                    if line.startswith('>'):
                        if current_id:
                            sequences[current_id] = ''.join(current_seq)
                        current_id = line[1:].split()[0]  # Take first word as ID
                        current_seq = []
                    elif line and current_id:
                        current_seq.append(line.upper())
                
                # Add last sequence
                if current_id and current_seq:
                    sequences[current_id] = ''.join(current_seq)
        
        except FileNotFoundError:
            print(f"Reference database file not found: {filepath}")
        
        return sequences
    
class MLFeatureExtractor:
    """Extract features from primer sets for machine learning"""
    
    @staticmethod
    def extract_primer_features(primer: Primer) -> np.ndarray:
        """Extract features from a single primer"""
        seq = primer.sequence.upper()
        
        features = [
            # Basic properties
            len(seq),
            primer.tm,
            primer.gc_content,
            primer.dg_3prime,
            primer.dg_hairpin,
            
            # Sequence composition
            seq.count('A') / len(seq),
            seq.count('T') / len(seq),
            seq.count('G') / len(seq),
            seq.count('C') / len(seq),
            
            # Dinucleotide frequencies
            sum(seq.count(dn) for dn in ['AA', 'AT', 'AC', 'AG']) / max(1, len(seq) - 1),
            sum(seq.count(dn) for dn in ['TA', 'TT', 'TC', 'TG']) / max(1, len(seq) - 1),
            sum(seq.count(dn) for dn in ['GA', 'GT', 'GC', 'GG']) / max(1, len(seq) - 1),
            sum(seq.count(dn) for dn in ['CA', 'CT', 'CC', 'CG']) / max(1, len(seq) - 1),
            
            # Structural features
            1 if seq.startswith(('G', 'C')) else 0,  # GC clamp start
            1 if seq.endswith(('G', 'C')) else 0,   # GC clamp end
            max(len(match.group()) for match in re.finditer(r'(.)\1+', seq)) if re.search(r'(.)\1+', seq) else 1,  # Max homopolymer
            
            # 3' end composition (last 5 bases)
            (seq[-5:].count('G') + seq[-5:].count('C')) / 5 if len(seq) >= 5 else 0,
        ]
        
        return np.array(features)
    
    @staticmethod
    def extract_primer_set_features(lamp_set: LAMPPrimerSet) -> np.ndarray:
        """Extract features from a complete primer set"""
        primers = [lamp_set.f3, lamp_set.b3, lamp_set.fip, lamp_set.bip]
        
        # Individual primer features
        primer_features = []
        for primer in primers:
            primer_features.extend(MLFeatureExtractor.extract_primer_features(primer))
        
        # Inter-primer features
        tms = [p.tm for p in primers]
        gcs = [p.gc_content for p in primers]
        dg_3primes = [p.dg_3prime for p in primers]
        
        inter_features = [
            # Tm statistics
            np.mean(tms),
            np.std(tms),
            max(tms) - min(tms),
            
            # GC statistics
            np.mean(gcs),
            np.std(gcs),
            max(gcs) - min(gcs),
            
            # ΔG statistics
            np.mean(dg_3primes),
            np.std(dg_3primes),
            
            # Primer lengths
            np.mean([len(p.sequence) for p in primers]),
            np.std([len(p.sequence) for p in primers]),
        ]
        
        return np.concatenate([primer_features, inter_features])

class MLPrimerPredictor:
    """Machine learning model for primer quality prediction"""
    
    def __init__(self):
        self.model = None
        self.scaler = None
        self.is_trained = False
        
    def prepare_training_data(self, primer_sets: List[LAMPPrimerSet], 
                            scores: List[float]) -> Tuple[np.ndarray, np.ndarray]:
        """Prepare training data from primer sets and their scores"""
        X = []
        y = np.array(scores)
        
        for lamp_set in primer_sets:
            features = MLFeatureExtractor.extract_primer_set_features(lamp_set)
            X.append(features)
        
        X = np.array(X)
        
        # Normalize scores to 0-1 range for better training
        if len(y) > 1:
            y_normalized = (y - np.min(y)) / (np.max(y) - np.min(y) + 1e-8)
        else:
            y_normalized = y
            
        return X, y_normalized
    
    def train_model(self, X: np.ndarray, y: np.ndarray, test_size: float = 0.2):
        """Train ML model on primer set features"""
        try:
            from sklearn.ensemble import RandomForestRegressor
            from sklearn.model_selection import train_test_split
            from sklearn.preprocessing import StandardScaler
            from sklearn.metrics import mean_squared_error, r2_score
            
            # Split data
            if len(X) > 10:  # Only split if we have enough data
                X_train, X_test, y_train, y_test = train_test_split(
                    X, y, test_size=test_size, random_state=42
                )
            else:
                X_train, X_test = X, X
                y_train, y_test = y, y
            
            # Scale features
            self.scaler = StandardScaler()
            X_train_scaled = self.scaler.fit_transform(X_train)
            X_test_scaled = self.scaler.transform(X_test)
            
            # Train Random Forest model
            self.model = RandomForestRegressor(
                n_estimators=100,
                max_depth=10,
                random_state=42,
                n_jobs=-1
            )
            
            self.model.fit(X_train_scaled, y_train)
            
            # Evaluate model
            y_pred = self.model.predict(X_test_scaled)
            mse = mean_squared_error(y_test, y_pred)
            r2 = r2_score(y_test, y_pred)
            
            self.is_trained = True
            
            return {
                'mse': mse,
                'r2': r2,
                'train_size': len(X_train),
                'test_size': len(X_test)
            }
            
        except ImportError:
            print("scikit-learn not available. Using rule-based scoring.")
            return None
    
    def predict_quality(self, lamp_set: LAMPPrimerSet) -> float:
        """Predict quality score for a primer set"""
        if not self.is_trained or self.model is None:
            return 0.0
        
        features = MLFeatureExtractor.extract_primer_set_features(lamp_set)
        features_scaled = self.scaler.transform(features.reshape(1, -1))
        
        return self.model.predict(features_scaled)[0]
    
    def get_feature_importance(self) -> Optional[np.ndarray]:
        """Get feature importance from trained model"""
        if self.is_trained and hasattr(self.model, 'feature_importances_'):
            return self.model.feature_importances_
        return None

@dataclass
class OffTargetMatch:
    """Represents an off-target match for a primer"""
    sequence_id: str
    position: int
    primer_sequence: str
    target_sequence: str
    tm: float
    mismatches: int
    gaps: int
    alignment_score: float

@dataclass
class BatchTarget:
    """Represents a target sequence for batch processing"""
    target_id: str
    sequence: str
    description: Optional[str] = None
    seq_type: Optional[SequenceType] = None
    custom_config: Optional[Dict] = None

@dataclass
class BatchResults:
    """Results from batch primer design"""
    target_results: Dict[str, DesignResults]
    batch_statistics: Dict[str, any]
    cross_target_analysis: Dict[str, any]
    batch_effects: Dict[str, any]
    
    def get_best_primers_per_target(self) -> Dict[str, LAMPPrimerSet]:
        """Get best primer set for each target"""
        best_primers = {}
        for target_id, results in self.target_results.items():
            if results.primer_sets:
                best_primers[target_id] = results.primer_sets[0]
        return best_primers
    
    def get_batch_summary(self) -> Dict:
        """Get summary statistics across all targets"""
        total_sets = sum(len(r.primer_sets) for r in self.target_results.values())
        successful_targets = sum(1 for r in self.target_results.values() if r.primer_sets)
        
        return {
            'total_targets': len(self.target_results),
            'successful_targets': successful_targets,
            'total_primer_sets': total_sets,
            'success_rate': successful_targets / len(self.target_results) if self.target_results else 0,
            'batch_effects': self.batch_effects,
            'cross_target_interactions': self.cross_target_analysis
        }
    """Results of specificity analysis"""
    primer: Primer
    off_targets: List[OffTargetMatch]
    is_specific: bool
    max_off_target_tm: float
    total_off_targets: int
    """Container for primer design results"""
    primer_sets: List[LAMPPrimerSet]
    total_candidates: int
    total_combinations: int
    filtering_stats: Dict[str, int]
    scores: List[float]
    
    def get_top_n(self, n: int = 3) -> List[LAMPPrimerSet]:
        """Get top N primer sets"""
        return self.primer_sets[:n]
    
    def get_all(self) -> List[LAMPPrimerSet]:
        """Get all valid primer sets"""
        return self.primer_sets
    
    def get_statistics(self) -> Dict:
        """Get design statistics"""
        return {
            'total_candidates': self.total_candidates,
            'total_combinations': self.total_combinations,
            'valid_sets': len(self.primer_sets),
            'filtering_stats': self.filtering_stats,
            'score_range': (min(self.scores) if self.scores else 0, 
                          max(self.scores) if self.scores else 0)
        }

class LAMPPrimerDesigner:
    """Main class for LAMP primer design with optional filtering control"""
    
    def __init__(self, config: PrimerConfig):
        self.config = config
        self.thermo = DNAThermodynamics()
        self.processor = SequenceProcessor()
        self.ml_predictor = MLPrimerPredictor() if config.use_ml_scoring else None
        self.specificity_checker = NonSpecificityChecker(config) if config.enable_specificity_check else None
        self.disable_filtering = False  # New flag for --full mode
        self.filtering_stats = {
            'total_candidates': 0,
            'tm_filtered': 0,
            'gc_filtered': 0,
            'dg_filtered': 0,
            'quality_filtered': 0,
            'distance_filtered': 0,
            'dimer_filtered': 0,
            'specificity_filtered': 0,
            'filtering_disabled': 0
        }

    def _validate_and_clean_sequence(self, sequence: str) -> str:
        """Permissive sequence validation and cleaning"""
        
        # Convert to uppercase
        sequence = sequence.upper().strip()
        
        # Remove any whitespace, numbers, or common FASTA artifacts
        sequence = re.sub(r'[^ATGCNRYSWKMBDHV]', '', sequence)
        
        # Basic validation - accept standard IUPAC nucleotide codes
        valid_chars = set('ATGCNRYSWKMBDHV')
        sequence_chars = set(sequence)
        
        if not sequence_chars.issubset(valid_chars):
            invalid_chars = sequence_chars - valid_chars
            print(f"⚠️  Warning: Removing invalid characters: {invalid_chars}")
            # Keep only valid characters
            sequence = ''.join(c for c in sequence if c in valid_chars)
        
        # Convert ambiguous bases to N for processing
        ambiguous_map = {
            'R': 'N', 'Y': 'N', 'S': 'N', 'W': 'N', 'K': 'N', 'M': 'N',
            'B': 'N', 'D': 'N', 'H': 'N', 'V': 'N'
        }
        
        cleaned_sequence = ''
        for char in sequence:
            if char in ambiguous_map:
                cleaned_sequence += ambiguous_map[char]
                if len(cleaned_sequence) <= 50:  # Only warn for first 50 characters
                    print(f"⚠️  Converting ambiguous base {char} to N at position {len(cleaned_sequence)}")
            else:
                cleaned_sequence += char
        
        return cleaned_sequence
    
    def find_primers_with_type_specific_thresholds(self, sequence: str, 
                                                 type_specific_designer: PrimerTypeSpecificDesigner) -> Tuple[Dict[LAMPPrimerType, List], Dict]:
        """Find primers using type-specific thresholds for each primer type."""
        
        # Analyze sequence and get adapted configs
        adapted_configs = type_specific_designer.analyze_and_adapt_all_primers(sequence)
        type_specific_designer.print_primer_type_summary()
        
        candidates = {primer_type: [] for primer_type in LAMPPrimerType if primer_type not in [LAMPPrimerType.FIP, LAMPPrimerType.BIP]}
        
        # Define search regions
        seq_len = len(sequence)
        regions = {
            LAMPPrimerType.F3: (0, seq_len // 2),
            LAMPPrimerType.F2: (seq_len // 6, 4 * seq_len // 6),
            LAMPPrimerType.F1C: (seq_len // 4, 3 * seq_len // 4),
            LAMPPrimerType.B1C: (seq_len // 4, 3 * seq_len // 4),
            LAMPPrimerType.B2: (2 * seq_len // 6, 5 * seq_len // 6),
            LAMPPrimerType.B3: (seq_len // 2, seq_len)
        }
        
        print(f"\n🔍 PRIMER-TYPE-SPECIFIC SEARCH:")
        print(f"{'='*50}")
        
        # Search for each primer type with its specific thresholds
        for primer_type, (start_region, end_region) in regions.items():
            if primer_type in [LAMPPrimerType.FIP, LAMPPrimerType.BIP]:
                continue
            
            config = adapted_configs[primer_type]
            primer_candidates = []
            
            print(f"\n🎯 Searching {primer_type.value}:")
            print(f"   Region: {start_region}-{end_region}")
            print(f"   Tm target: {config['tm_optimal']:.1f}°C (range: {config['tm_range'][0]:.1f}-{config['tm_range'][1]:.1f}°C)")
            print(f"   Length: {config['length_range'][0]}-{config['length_range'][1]} bp")
            
            valid_count = 0
            total_tested = 0
            
            # Search in the region
            for start in range(start_region, min(end_region, seq_len - config['length_range'][0]), 2):
                for length in range(config['length_range'][0], min(config['length_range'][1] + 1, seq_len - start + 1)):
                    primer_seq = sequence[start:start + length]
                    total_tested += 1
                    
                    # Validate against type-specific config
                    is_valid, properties = self.validate_primer_against_type_config(
                        primer_seq, primer_type, config
                    )
                    
                    if is_valid:
                        valid_count += 1
                        primer_candidates.append((properties, properties['quality_score']))
            
            # Sort by quality score and keep best candidates
            primer_candidates.sort(key=lambda x: x[1], reverse=True)
            candidates[primer_type] = [props for props, score in primer_candidates[:30]]
            
            print(f"   Found: {valid_count}/{total_tested} valid candidates")
            if candidates[primer_type]:
                best = candidates[primer_type][0]
                print(f"   Best: {best['sequence'][:20]}... (Tm={best['tm']:.1f}°C, Score={best['quality_score']:.3f})")
            
        return candidates, adapted_configs
    
    def create_balanced_lamp_sets_with_type_awareness(self, candidates: Dict, 
                                                    adapted_configs: Dict) -> List:
        """Create LAMP primer sets considering type-specific balances and constraints."""
        
        lamp_sets = []
        max_combinations = 50000
        combinations_tested = 0
        
        print(f"\n🔄 CREATING TYPE-AWARE LAMP COMBINATIONS:")
        print(f"{'='*50}")
        
        # Get primer type configs for validation
        f3_config = adapted_configs[LAMPPrimerType.F3]
        b3_config = adapted_configs[LAMPPrimerType.B3]
        fip_config = adapted_configs[LAMPPrimerType.FIP]
        bip_config = adapted_configs[LAMPPrimerType.BIP]
        
        # Limits for each primer type
        limits = {pt: min(15, len(candidates.get(pt, []))) for pt in candidates.keys()}
        
        for f3 in candidates.get(LAMPPrimerType.F3, [])[:limits.get(LAMPPrimerType.F3, 0)]:
            for b3 in candidates.get(LAMPPrimerType.B3, [])[:limits.get(LAMPPrimerType.B3, 0)]:
                
                # Check outer primer balance first
                outer_tm_diff = abs(f3['tm'] - b3['tm'])
                if outer_tm_diff > 5.0:
                    continue
                    
                for f1c in candidates.get(LAMPPrimerType.F1C, [])[:limits.get(LAMPPrimerType.F1C, 0)]:
                    for f2 in candidates.get(LAMPPrimerType.F2, [])[:limits.get(LAMPPrimerType.F2, 0)]:
                        for b1c in candidates.get(LAMPPrimerType.B1C, [])[:limits.get(LAMPPrimerType.B1C, 0)]:
                            for b2 in candidates.get(LAMPPrimerType.B2, [])[:limits.get(LAMPPrimerType.B2, 0)]:
                                
                                combinations_tested += 1
                                if combinations_tested > max_combinations:
                                    break
                                
                                # Create composite primers
                                fip_seq = self.processor.reverse_complement(f1c['sequence']) + f2['sequence']
                                bip_seq = self.processor.reverse_complement(b1c['sequence']) + b2['sequence']
                                
                                # Calculate composite primer properties
                                fip_tm = (f1c['tm'] + f2['tm']) / 2
                                bip_tm = (b1c['tm'] + b2['tm']) / 2
                                fip_gc = self.thermo.calculate_gc_content(fip_seq)
                                bip_gc = self.thermo.calculate_gc_content(bip_seq)
                                
                                # Validate composite primers
                                fip_valid = (fip_config['tm_range'][0] <= fip_tm <= fip_config['tm_range'][1] and
                                        fip_config['gc_range'][0] <= fip_gc <= fip_config['gc_range'][1] and
                                        fip_config['length_range'][0] <= len(fip_seq) <= fip_config['length_range'][1])
                                
                                bip_valid = (bip_config['tm_range'][0] <= bip_tm <= bip_config['tm_range'][1] and
                                        bip_config['gc_range'][0] <= bip_gc <= bip_config['gc_range'][1] and
                                        bip_config['length_range'][0] <= len(bip_seq) <= bip_config['length_range'][1])
                                
                                if not (fip_valid and bip_valid):
                                    continue
                                
                                # Check hierarchical Tm relationships
                                all_tms = [f3['tm'], b3['tm'], fip_tm, bip_tm]
                                tm_range = max(all_tms) - min(all_tms)
                                
                                # Type-aware Tm balance
                                outer_avg = (f3['tm'] + b3['tm']) / 2
                                inner_avg = (fip_tm + bip_tm) / 2
                                
                                if outer_avg < inner_avg - 2.0:
                                    continue
                                
                                if tm_range <= 8.0:
                                    lamp_set = {
                                        'f3': f3,
                                        'b3': b3, 
                                        'f1c': f1c,
                                        'f2': f2,
                                        'b1c': b1c,
                                        'b2': b2,
                                        'fip_seq': fip_seq,
                                        'bip_seq': bip_seq,
                                        'fip_tm': fip_tm,
                                        'bip_tm': bip_tm,
                                        'fip_gc': fip_gc,
                                        'bip_gc': bip_gc,
                                        'tm_range': tm_range,
                                        'tm_hierarchy_score': outer_avg - inner_avg,
                                        'overall_quality': np.mean([f3['quality_score'], b3['quality_score'], 
                                                                f1c['quality_score'], f2['quality_score'],
                                                                b1c['quality_score'], b2['quality_score']])
                                    }
                                    
                                    lamp_sets.append(lamp_set)
                                
                                if combinations_tested >= max_combinations:
                                    break
                            if combinations_tested >= max_combinations:
                                break
                        if combinations_tested >= max_combinations:
                            break
                    if combinations_tested >= max_combinations:
                        break
                if combinations_tested >= max_combinations:
                    break
            if combinations_tested >= max_combinations:
                break
        
        print(f"   Tested: {combinations_tested:,} combinations")
        print(f"   Valid sets: {len(lamp_sets)}")
        
        # Sort by combined quality and hierarchy score
        lamp_sets.sort(key=lambda x: (x['overall_quality'] + x['tm_hierarchy_score'] * 0.1), reverse=True)
        
        if lamp_sets:
            best = lamp_sets[0]
            print(f"   Best Tm range: {best['tm_range']:.1f}°C")
            print(f"   Best hierarchy score: {best['tm_hierarchy_score']:.1f}°C")
            print(f"   Best overall quality: {best['overall_quality']:.3f}")
        
        return lamp_sets
    
    def design_primers(self, target_sequence: str, return_format: ResultFormat = ResultFormat.TOP_N, 
                          top_n: int = 3, reference_sequences: Optional[Dict[str, str]] = None) -> Union[List[LAMPPrimerSet], DesignResults]:
        """Enhanced LAMP primer design with flexible Tm-balanced approach"""
        try:
            # Validate input length
            if not target_sequence or len(target_sequence) < 100:
                raise ValueError("Target sequence too short for LAMP design (minimum 100bp)")
            
            # FIXED: More permissive sequence validation
            processed_seq = self._validate_and_clean_sequence(target_sequence)
            
            if len(processed_seq) < 100:
                raise ValueError(f"Cleaned sequence too short: {len(processed_seq)}bp (minimum 100bp)")
            
            print(f"✅ Sequence validation passed: {len(processed_seq)}bp cleaned sequence")
            
            # Process sequence based on type
            if self.config.seq_type == SequenceType.BISULFITE:
                processed_seq = self.processor.process_bisulfite(processed_seq)
            
            # Load reference database for specificity checking (only if filtering enabled)
            if not self.disable_filtering and self.specificity_checker and reference_sequences:
                self.specificity_checker.load_reference_database(sequences=reference_sequences)
            elif not self.disable_filtering and self.specificity_checker and self.config.specificity_database:
                self.specificity_checker.load_reference_database(self.config.specificity_database)
            
            # Analyze input sequence
            self.analyze_input_sequence(processed_seq)

            if self.disable_filtering:
                print("Note: Quality filtering disabled")
            else:
                self.debug_sequence_analysis(processed_seq)
            
            # Find potential primer regions using FLEXIBLE approach
            primer_candidates = self._find_primer_candidates_flexible(processed_seq)
            
            # Generate LAMP primer combinations with Tm balance
            lamp_sets = self._generate_lamp_combinations_balanced(primer_candidates, processed_seq)
            
            if len(lamp_sets) == 0:
                print("⚠️  Standard generation failed, trying emergency mode...")
                lamp_sets = self.emergency_primer_generation(processed_seq, primer_candidates)
            elif len(lamp_sets) == 0:
                print(f"\n❌ No primer sets found after all attempts")
                self.analyze_why_no_primers(processed_seq, primer_candidates)
                print(f"\n💡 SUGGESTIONS:")
                print(f"   1. Try with --full flag to disable all filtering")
                print(f"   2. Relax Tm constraints (try -t 45 -T 85)")
                print(f"   3. Relax GC constraints (try -c 15 -C 85)")
                print(f"   4. Your sequence is {len(processed_seq)}bp - this should work!")

            # Filter and rank primer sets (skip if filtering disabled)
            if self.disable_filtering:
                print(f"⚠️  Filtering disabled: Returning {len(lamp_sets)} raw primer combinations")
                valid_sets = lamp_sets[:100]  # Limit for performance
                self.filtering_stats['filtering_disabled'] = len(lamp_sets)
                specificity_results = None
            else:
                valid_sets = self._filter_primer_sets_balanced(lamp_sets)
                
                # Check specificity if enabled
                specificity_results = None
                if self.specificity_checker and valid_sets:
                    print(f"Checking specificity for {len(valid_sets)} primer sets...")
                    valid_sets, specificity_results = self._filter_by_specificity(valid_sets)
                    print(f"After specificity filtering: {len(valid_sets)} primer sets remain")
            
            # Train ML model if enabled and we have enough data
            if self.config.use_ml_scoring and len(valid_sets) >= 50:
                self._train_ml_model(valid_sets)
            
            # Score and sort primer sets
            scored_sets = self._score_and_sort_primer_sets(valid_sets)
            
            # Prepare results based on format
            if return_format == ResultFormat.TOP_N:
                return scored_sets[:top_n]
            elif return_format == ResultFormat.ALL:
                return scored_sets
            else:  # DETAILED
                scores = [self._score_primer_set(lamp_set) for lamp_set in scored_sets]
                return DesignResults(
                    primer_sets=scored_sets,
                    total_candidates=self.filtering_stats['total_candidates'],
                    total_combinations=len(lamp_sets),
                    filtering_stats=self.filtering_stats.copy(),
                    scores=scores,
                    specificity_results=specificity_results
                )
            
        except Exception as e:
            print(f"❌ Error in primer design: {e}")
            print(f"   Sequence length: {len(target_sequence) if target_sequence else 0}")
            print(f"   First 100 chars: {target_sequence[:100] if target_sequence else 'None'}")
            
            # Return empty results instead of crashing
            if return_format == ResultFormat.DETAILED:
                return DesignResults(
                    primer_sets=[],
                    total_candidates=0,
                    total_combinations=0,
                    filtering_stats=self.filtering_stats.copy(),
                    scores=[]
                )
            else:
                return []
            
    def debug_filtering(self, primer_seq: str, tm: float, gc: float, 
                   dg_3prime: float, dg_hairpin: float) -> str:
        """Debug helper to see why primers are being filtered"""
        reasons = []
        
        if not (self.config.min_tm <= tm <= self.config.max_tm):
            reasons.append(f"Tm {tm:.1f} not in [{self.config.min_tm}, {self.config.max_tm}]")
        
        if not (self.config.min_gc <= gc <= self.config.max_gc):
            reasons.append(f"GC {gc:.1f}% not in [{self.config.min_gc}, {self.config.max_gc}]")
        
        if not (dg_3prime >= self.config.max_3prime_dg):
            reasons.append(f"3' ΔG {dg_3prime:.2f} < {self.config.max_3prime_dg}")
        
        if not (dg_hairpin >= self.config.max_hairpin_dg):
            reasons.append(f"Hairpin ΔG {dg_hairpin:.2f} < {self.config.max_hairpin_dg}")
        
        return "; ".join(reasons) if reasons else "PASSED"
    
    def debug_candidate_search(self, sequence: str, primer_type: LAMPPrimerType, config: Dict):
        """Debug why type-specific search finds no candidates"""
        
        seq_len = len(sequence)
        
        print(f"\n🔍 DEBUGGING {primer_type.value}:")
        print(f"   Required Tm: {config['tm_range'][0]:.1f}-{config['tm_range'][1]:.1f}°C")
        print(f"   Required length: {config['length_range'][0]}-{config['length_range'][1]}bp")
        print(f"   Required GC: {config['gc_range'][0]:.1f}-{config['gc_range'][1]:.1f}%")
        
        # Test samples across the sequence
        sample_positions = [50, seq_len//4, seq_len//2, 3*seq_len//4, seq_len-50]
        
        found_valid = False
        for pos in sample_positions:
            if pos + 20 <= seq_len:
                sample_seq = sequence[pos:pos+20]
                tm = DNAThermodynamics.calculate_tm(sample_seq)
                gc = DNAThermodynamics.calculate_gc_content(sample_seq)
                
                tm_ok = config['tm_range'][0] <= tm <= config['tm_range'][1]
                gc_ok = config['gc_range'][0] <= gc <= config['gc_range'][1]
                
                status = "✅" if (tm_ok and gc_ok) else "❌"
                print(f"   Sample at {pos}: {sample_seq[:15]}... Tm={tm:.1f}°C, GC={gc:.1f}% {status}")
                
                if tm_ok and gc_ok:
                    found_valid = True
        
        if not found_valid:
            print(f"   ❌ No samples pass - constraints may be too strict!")
            # Suggest relaxed constraints
            all_tms = []
            all_gcs = []
            for pos in range(0, seq_len-20, 20):
                if pos + 20 <= seq_len:
                    sample_seq = sequence[pos:pos+20]
                    all_tms.append(DNAThermodynamics.calculate_tm(sample_seq))
                    all_gcs.append(DNAThermodynamics.calculate_gc_content(sample_seq))
            
            if all_tms and all_gcs:
                print(f"   💡 Actual Tm range in sequence: {min(all_tms):.1f}-{max(all_tms):.1f}°C")
                print(f"   💡 Actual GC range in sequence: {min(all_gcs):.1f}-{max(all_gcs):.1f}%")
    
    def debug_sequence_analysis(self, sequence: str):
        """Detailed analysis when no primers are found"""
        
        print(f"\n🔬 DETAILED SEQUENCE ANALYSIS:")
        print(f"{'='*60}")
        
        # Check for problematic regions
        print(f"\n📊 Sequence Features:")
        
        # Homopolymer runs
        for base in ['A', 'T', 'G', 'C']:
            matches = re.finditer(f'{base}{{4,}}', sequence)
            runs = [(m.start(), m.end(), m.end() - m.start()) for m in matches]
            if runs:
                print(f"   {base}-runs: {len(runs)} found")
                for start, end, length in runs[:3]:  # Show first 3
                    print(f"      Position {start}-{end}: {base*length}")
        
        # GC distribution
        window_size = 100
        gc_windows = []
        for i in range(0, len(sequence) - window_size, window_size // 2):
            window = sequence[i:i + window_size]
            gc = self.thermo.calculate_gc_content(window)
            gc_windows.append((i, gc))
        
        print(f"\n   GC Content Distribution (100bp windows):")
        for pos, gc in gc_windows[::4]:  # Show every 4th window
            print(f"      Position {pos}: {gc:.1f}%")
        
        # Test primer extraction at specific positions
        print(f"\n   Test Primer Extraction:")
        test_positions = [50, 200, 400, 600, 800, 1000]
        
        for pos in test_positions:
            if pos + 20 <= len(sequence):
                test_seq = sequence[pos:pos + 20]
                tm = self.thermo.calculate_tm(test_seq)
                gc = self.thermo.calculate_gc_content(test_seq)
                quality = self._check_primer_quality(test_seq)
                quality_flex = self._check_primer_quality_flexible(test_seq)
                
                print(f"      Pos {pos}: Tm={tm:.1f}°C, GC={gc:.1f}%, Quality={quality}, Flex={quality_flex}")
                if not quality and quality_flex:
                    print(f"         → Would pass with flexible quality check!")

    def _find_primer_candidates_flexible(self, sequence: str) -> Dict[str, List[Primer]]:
        """Find primers with progressive relaxation of constraints"""
    
        candidates = {
            'F3': [], 'B3': [], 'F1c': [], 'F2': [], 'B1c': [], 'B2': []
        }
        
        seq_len = len(sequence)
        total_candidates = 0
        
        # Define relaxation levels
        relaxation_levels = [
            {'name': 'Standard', 'tm_relax': 0, 'gc_relax': 0, 'dg_relax': 0},
            {'name': 'Relaxed', 'tm_relax': 5, 'gc_relax': 10, 'dg_relax': 0.5},
            {'name': 'Very Relaxed', 'tm_relax': 10, 'gc_relax': 20, 'dg_relax': 1.0},
            {'name': 'Maximum', 'tm_relax': 15, 'gc_relax': 30, 'dg_relax': 2.0}
        ]
        
        print(f"\n🔍 Progressive primer search with relaxation:")
        
        # Try each relaxation level
        for level in relaxation_levels:
            level_found = False
            
            print(f"\n   Level: {level['name']}")
            print(f"   Tm: {self.config.min_tm - level['tm_relax']:.1f} - {self.config.max_tm + level['tm_relax']:.1f}°C")
            print(f"   GC: {max(5, self.config.min_gc - level['gc_relax']):.1f} - {min(95, self.config.max_gc + level['gc_relax']):.1f}%")
            
            # Define flexible regions
            flexible_regions = {
                'F3': (0, seq_len // 2),
                'F2': (seq_len // 6, 4 * seq_len // 6),
                'F1c': (seq_len // 4, 3 * seq_len // 4),
                'B1c': (seq_len // 4, 3 * seq_len // 4),
                'B2': (2 * seq_len // 6, 5 * seq_len // 6),
                'B3': (seq_len // 2, seq_len)
            }
            
            for primer_type, (start_region, end_region) in flexible_regions.items():
                if len(candidates[primer_type]) >= 20:  # Already have enough
                    continue
                    
                region_candidates = []
                step_size = 2  # Use step size 2
                
                for start in range(start_region, min(end_region, seq_len - self.config.min_length), step_size):
                    for length in range(self.config.min_length, 
                                    min(self.config.max_length + 1, seq_len - start + 1)):
                        end = start + length
                        primer_seq = sequence[start:end]
                        total_candidates += 1
                        
                        # Calculate properties
                        tm = self.thermo.calculate_tm(primer_seq, self.config.na_conc, self.config.mg_conc)
                        gc = self.thermo.calculate_gc_content(primer_seq)
                        dg_3prime = self.thermo.calculate_3prime_stability(primer_seq)
                        dg_hairpin = self.thermo.calculate_hairpin_dg(primer_seq)
                        
                        # Apply relaxed constraints
                        tm_ok = (self.config.min_tm - level['tm_relax'] <= tm <= 
                                self.config.max_tm + level['tm_relax'])
                        gc_ok = (max(5, self.config.min_gc - level['gc_relax']) <= gc <= 
                                min(95, self.config.max_gc + level['gc_relax']))
                        dg_3prime_ok = dg_3prime >= (self.config.max_3prime_dg - level['dg_relax'])
                        dg_hairpin_ok = dg_hairpin >= (self.config.max_hairpin_dg - level['dg_relax'])
                        
                        # Use flexible quality check
                        quality_ok = self._check_primer_quality_flexible(primer_seq)
                        
                        if tm_ok and gc_ok and dg_3prime_ok and dg_hairpin_ok and quality_ok:
                            primer = Primer(primer_seq, start, end, tm, gc, primer_type,
                                        dg_3prime, dg_hairpin)
                            region_candidates.append(primer)
                            level_found = True
                
                # Add new candidates
                candidates[primer_type].extend(region_candidates)
                if region_candidates:
                    print(f"      {primer_type}: +{len(region_candidates)} candidates (total: {len(candidates[primer_type])})")
            
            # Stop if we have enough candidates for all types
            if all(len(candidates[pt]) >= 10 for pt in candidates.keys()):
                print(f"   ✅ Sufficient candidates found at {level['name']} level")
                break
        
        self.filtering_stats['total_candidates'] = total_candidates
        
        # Final summary
        print(f"\n📊 Progressive Search Summary:")
        for primer_type, primer_list in candidates.items():
            if primer_list:
                tms = [p.tm for p in primer_list]
                gcs = [p.gc_content for p in primer_list]
                print(f"   {primer_type}: {len(primer_list)} candidates")
                print(f"      Tm range: {min(tms):.1f}-{max(tms):.1f}°C")
                print(f"      GC range: {min(gcs):.1f}-{max(gcs):.1f}%")
            else:
                print(f"   {primer_type}: 0 candidates ⚠️")
        
        return candidates
    
    def analyze_input_sequence(self, sequence: str):
        """Analyze input sequence to understand why no primers are found"""
        print(f"\n🔍 INPUT SEQUENCE ANALYSIS:")
        print(f"   Length: {len(sequence)} bp")
        print(f"   GC Content: {DNAThermodynamics.calculate_gc_content(sequence):.1f}%")
        print(f"   First 50 bp: {sequence[:50]}")
        print(f"   Last 50 bp: {sequence[-50:]}")
        
        # Check if sequence is too short
        min_total_length = self.config.max_length * 6 + 200  # Rough estimate
        if len(sequence) < min_total_length:
            print(f"   ⚠️  WARNING: Sequence might be too short for LAMP primer design")
            print(f"   ⚠️  Recommended minimum: ~{min_total_length} bp")
        
        # Sample a few potential primer sequences
        print(f"\n   Sample potential primers:")
        for i, region_name in enumerate(['5\' region', 'middle', '3\' region']):
            start_pos = i * len(sequence) // 3
            for length in [18, 25, 30]:
                if start_pos + length <= len(sequence):
                    sample_seq = sequence[start_pos:start_pos + length]
                    tm = DNAThermodynamics.calculate_tm(sample_seq)
                    gc = DNAThermodynamics.calculate_gc_content(sample_seq)
                    print(f"     {region_name} ({length}bp): {sample_seq[:20]}... Tm={tm:.1f}°C, GC={gc:.1f}%")

    def _check_flexible_positioning(self, f3, b3, f1c, f2, b1c, b2) -> bool:
        """Check flexible positioning constraints with more tolerance"""
        
        # Basic requirement: maintain general 5' to 3' order with flexibility
        # 1. F3 should be before B3
        if f3.start >= b3.start:
            return False
        
        # 2. Minimum amplicon size (very permissive)
        if b3.end - f3.start < 100:  # Absolute minimum for LAMP
            return False
        
        # 3. Maximum amplicon size
        if b3.end - f3.start > 1000:  # Too large for efficient LAMP
            return False
        
        # That's it! No other positioning constraints
        return True
    
    def _check_permissive_positioning(self, f3, b3, f1c, f2, b1c, b2, sequence) -> bool:
        """Very permissive positioning check focusing on LAMP functionality"""
        
        # Minimum requirements for LAMP to work
        min_amplicon_size = 80  # Minimum total amplicon
        max_amplicon_size = 1000  # Maximum for efficiency
        min_loop_size = 10       # Minimum loop regions
        
        # Get all primer positions
        all_starts = [f3.start, f2.start, f1c.start, b1c.start, b2.start, b3.start]
        all_ends = [f3.end, f2.end, f1c.end, b1c.end, b2.end, b3.end]
        
        # Calculate approximate amplicon
        amplicon_start = min(all_starts)
        amplicon_end = max(all_ends)
        amplicon_size = amplicon_end - amplicon_start
        
        print(f"      Position check: F3({f3.start}-{f3.end}), B3({b3.start}-{b3.end}), Amplicon: {amplicon_size}bp")
        
        # Check amplicon size (VERY permissive)
        if not (min_amplicon_size <= amplicon_size <= max_amplicon_size):
            print(f"      ❌ Amplicon size {amplicon_size} not in range [{min_amplicon_size}, {max_amplicon_size}]")
            return False
        
        # Very basic ordering - just ensure F3 is generally before B3 
        if f3.start >= b3.start:
            print(f"      ❌ F3 start ({f3.start}) >= B3 start ({b3.start})")
            return False
        
        # Ensure some separation between outer primers (VERY permissive)
        f3_b3_gap = b3.start - f3.end
        if f3_b3_gap < 30:  # Reduced from 50
            print(f"      ❌ F3-B3 gap {f3_b3_gap} < 30bp")
            return False
        
        print(f"      ✅ Position check passed")
        return True
    
    def _check_permissive_dimers(self, primers: List[Primer]) -> bool:
        """Permissive dimer checking - only reject very strong dimers"""
        strong_dimer_threshold = -8.0  # Only reject very strong dimers (was -3.0)
        
        for i, primer1 in enumerate(primers):
            for j, primer2 in enumerate(primers[i+1:], i+1):
                dimer_dg = self.thermo.calculate_dimer_dg(primer1.sequence, primer2.sequence)
                if dimer_dg < strong_dimer_threshold:
                    return False  # Only reject very strong dimers
        return True
    
    def _check_permissive_dimers_threshold(self, primers: List[Primer], threshold: float) -> bool:
        """Check dimers with specified threshold"""
        for i, primer1 in enumerate(primers):
            for j, primer2 in enumerate(primers[i+1:], i+1):
                dimer_dg = self.thermo.calculate_dimer_dg(primer1.sequence, primer2.sequence)
                if dimer_dg < threshold:
                    return False
        return True

    def _is_duplicate_set(self, new_set: LAMPPrimerSet, existing_sets: List[LAMPPrimerSet]) -> bool:
        """Check if this primer set is essentially a duplicate"""
        for existing in existing_sets:
            if (new_set.f3.sequence == existing.f3.sequence and
                new_set.b3.sequence == existing.b3.sequence and
                new_set.fip.sequence == existing.fip.sequence and
                new_set.bip.sequence == existing.bip.sequence):
                return True
        return False
    
    def _score_primer_positioning(self, lamp_set: LAMPPrimerSet, sequence: str) -> float:
        """Score primer positioning for LAMP efficiency"""
        
        # Calculate key distances
        f3_to_f2_dist = abs(lamp_set.fip.start - lamp_set.f3.end)  # Approximate F3 to F2 distance
        b3_to_b2_dist = abs(lamp_set.b3.start - lamp_set.bip.end)  # Approximate B3 to B2 distance
        
        # Ideal distances for LAMP (literature values)
        ideal_f3_f2 = 60  # 40-80bp typical
        ideal_b3_b2 = 60  # 40-80bp typical
        
        # Distance scoring (closer to ideal = higher score)
        f3_f2_score = 1.0 / (1.0 + abs(f3_to_f2_dist - ideal_f3_f2) / 20.0)
        b3_b2_score = 1.0 / (1.0 + abs(b3_to_b2_dist - ideal_b3_b2) / 20.0)
        
        # Amplicon size scoring
        amplicon_size = lamp_set.b3.end - lamp_set.f3.start
        ideal_amplicon = 300  # 200-400bp typical
        amplicon_score = 1.0 / (1.0 + abs(amplicon_size - ideal_amplicon) / 100.0)
        
        # Weighted positioning score
        positioning_score = (0.3 * f3_f2_score + 0.3 * b3_b2_score + 0.4 * amplicon_score)
        
        return positioning_score
    
    def _filter_primer_sets_balanced(self, lamp_sets: List[LAMPPrimerSet]) -> List[LAMPPrimerSet]:
        """Permissive filtering that only removes clearly problematic sets"""
        if self.disable_filtering:
            return lamp_sets
        
        valid_sets = []
        
        # VERY permissive Tm range - much larger than before
        max_allowed_tm_range = 12.0  # Increased from 3.0°C to 12.0°C
        
        print(f"\n🎯 Applying permissive Tm filter (max range: {max_allowed_tm_range}°C)...")
        
        for lamp_set in lamp_sets:
            tms = [lamp_set.f3.tm, lamp_set.b3.tm, lamp_set.fip.tm, lamp_set.bip.tm]
            tm_range = max(tms) - min(tms)
            
            # Only filter out sets with extremely bad Tm balance
            if tm_range <= max_allowed_tm_range:
                valid_sets.append(lamp_set)
        
        print(f"   Permissive Tm filter: {len(valid_sets)}/{len(lamp_sets)} sets passed")
        
        # If still no sets, try even more permissive
        if len(valid_sets) == 0 and len(lamp_sets) > 0:
            print(f"   No sets passed - using most balanced available sets...")
            # Sort by Tm balance and take best ones
            lamp_sets_with_balance = []
            for lamp_set in lamp_sets:
                tms = [lamp_set.f3.tm, lamp_set.b3.tm, lamp_set.fip.tm, lamp_set.bip.tm]
                tm_range = max(tms) - min(tms)
                lamp_sets_with_balance.append((lamp_set, tm_range))
            
            # Sort by Tm range (smaller is better)
            lamp_sets_with_balance.sort(key=lambda x: x[1])
            
            # Take the best 10 sets even if they don't meet strict criteria
            valid_sets = [lamp_set for lamp_set, tm_range in lamp_sets_with_balance[:30]]
            print(f"   Rescued {len(valid_sets)} best available sets")
            
            if valid_sets:
                best_range = lamp_sets_with_balance[0][1]
                print(f"   Best rescued Tm range: {best_range:.1f}°C")
        
        return valid_sets

    def analyze_why_no_primers(self, sequence: str, candidates: Dict[str, List[Primer]]):
        """Detailed analysis of why no primers were found"""
        print(f"\n🔍 DEBUGGING: Why no valid primer combinations?")
        
        seq_len = len(sequence)
        print(f"   Sequence length: {seq_len} bp")
        
        # Check candidate distribution
        print(f"\n   📊 Candidate positions:")
        for primer_type, primer_list in candidates.items():
            if primer_list:
                starts = [p.start for p in primer_list[:5]]  # First 5 positions
                tms = [p.tm for p in primer_list[:5]]
                print(f"   {primer_type}: positions {starts}, Tm range {min(tms):.1f}-{max(tms):.1f}°C")
            else:
                print(f"   {primer_type}: No candidates!")
        
        # Test a simple combination manually
        if all(len(candidates[pt]) > 0 for pt in ['F3', 'B3', 'F1c', 'F2', 'B1c', 'B2']):
            print(f"\n   🧪 Testing simple combination:")
            f3 = candidates['F3'][0]
            b3 = candidates['B3'][0]
            f1c = candidates['F1c'][0]
            f2 = candidates['F2'][0]
            b1c = candidates['B1c'][0]
            b2 = candidates['B2'][0]
            
            print(f"   F3: pos {f3.start}-{f3.end}, Tm={f3.tm:.1f}°C")
            print(f"   B3: pos {b3.start}-{b3.end}, Tm={b3.tm:.1f}°C")
            print(f"   F1c: pos {f1c.start}-{f1c.end}, Tm={f1c.tm:.1f}°C")
            print(f"   F2: pos {f2.start}-{f2.end}, Tm={f2.tm:.1f}°C")
            
            # Check positioning
            positioning_ok = self._check_permissive_positioning(f3, b3, f1c, f2, b1c, b2, sequence)
            print(f"   Positioning OK: {positioning_ok}")
            
            if not positioning_ok:
                print(f"   ❌ Position analysis:")
                print(f"      F3 end: {f3.end}, B3 start: {b3.start}")
                print(f"      F3-B3 distance: {b3.start - f3.end}")
                print(f"      Total span: {b3.end - f3.start} bp")
            
            # Check Tm balance
            all_tms = [f3.tm, b3.tm, (f1c.tm + f2.tm)/2, (b1c.tm + b2.tm)/2]
            tm_range = max(all_tms) - min(all_tms)
            print(f"   Tm range: {tm_range:.1f}°C (F3:{f3.tm:.1f}, B3:{b3.tm:.1f}, FIP:{(f1c.tm + f2.tm)/2:.1f}, BIP:{(b1c.tm + b2.tm)/2:.1f})")
            
            if tm_range > 15:
                print(f"   ❌ Tm range too large - sequence may need different parameters")

    def _generate_lamp_combinations_balanced(self, candidates: Dict[str, List[Primer]], 
                                         sequence: str) -> List[LAMPPrimerSet]:
        """Generate LAMP combinations with progressive constraint relaxation"""
        # Analyze sequence composition
        gc_content = self.thermo.calculate_gc_content(sequence)
        print(f"\n🧬 Sequence GC: {gc_content:.1f}% - using ULTRA-permissive mode")
        
        lamp_sets = []
        combinations_tested = 0
        max_combinations = 1000000  # Increased limit
        
        # Use more candidates
        limit_per_type = min(30, max(10, min(len(candidates[pt]) for pt in candidates.keys())))
        
        print(f"\n🔄 Generating combinations (up to {limit_per_type} candidates per type)...")
        print(f"   Using ULTRA-permissive mode due to challenging sequence")
        
        # Start with VERY permissive Tm tolerance
        tm_tolerances = [20.0, 25.0, 30.0, 35.0, 40.0]  # Much more permissive!
        
        for tm_tolerance in tm_tolerances:
            print(f"\n   Trying ULTRA-permissive Tm tolerance: ±{tm_tolerance}°C...")
            initial_count = len(lamp_sets)
            
            # Try all combinations
            for f3_idx, f3 in enumerate(candidates['F3'][:limit_per_type]):
                for b3_idx, b3 in enumerate(candidates['B3'][:limit_per_type]):
                    
                    # Very basic outer primer check
                    if abs(f3.tm - b3.tm) > tm_tolerance * 1.5:  # Even more permissive for outer primers
                        continue
                        
                    for f1c in candidates['F1c'][:limit_per_type]:
                        for f2 in candidates['F2'][:limit_per_type]:
                            for b1c in candidates['B1c'][:limit_per_type]:
                                for b2 in candidates['B2'][:limit_per_type]:
                                    combinations_tested += 1
                                    
                                    if combinations_tested % 100000 == 0:
                                        print(f"      Tested {combinations_tested:,} combinations...")
                                    
                                    if combinations_tested > max_combinations:
                                        break
                                    
                                    # Ultra-simple positioning check
                                    if not self._check_flexible_positioning(f3, b3, f1c, f2, b1c, b2):
                                        continue
                                    
                                    # Create composite primers
                                    fip_seq = self.processor.reverse_complement(f1c.sequence) + f2.sequence
                                    bip_seq = self.processor.reverse_complement(b1c.sequence) + b2.sequence
                                    
                                    # Calculate properties with error handling
                                    try:
                                        fip_tm = self.thermo.calculate_tm(fip_seq)
                                        bip_tm = self.thermo.calculate_tm(bip_seq)
                                        fip_gc = self.thermo.calculate_gc_content(fip_seq)
                                        bip_gc = self.thermo.calculate_gc_content(bip_seq)
                                        
                                        # Create primers
                                        fip = Primer(fip_seq, f1c.start, f2.end, fip_tm, fip_gc, 'FIP',
                                                self.thermo.calculate_3prime_stability(fip_seq),
                                                self.thermo.calculate_hairpin_dg(fip_seq))
                                        
                                        bip = Primer(bip_seq, b1c.start, b2.end, bip_tm, bip_gc, 'BIP',
                                                self.thermo.calculate_3prime_stability(bip_seq),
                                                self.thermo.calculate_hairpin_dg(bip_seq))
                                    except:
                                        continue  # Skip if calculation fails
                                    
                                    # Check Tm balance with ULTRA-permissive tolerance
                                    all_tms = [f3.tm, b3.tm, fip.tm, bip.tm]
                                    tm_range = max(all_tms) - min(all_tms)
                                    
                                    if tm_range <= tm_tolerance:
                                        # ULTRA-permissive dimer check
                                        lamp_set = LAMPPrimerSet(f3, b3, fip, bip)
                                        
                                        # Basic duplicate check
                                        is_dup = False
                                        for existing in lamp_sets:
                                            if (lamp_set.f3.sequence == existing.f3.sequence and
                                                lamp_set.b3.sequence == existing.b3.sequence):
                                                is_dup = True
                                                break
                                        
                                        if not is_dup:
                                            lamp_sets.append(lamp_set)
                                            print(f"      ✅ Found set #{len(lamp_sets)}: Tm range {tm_range:.1f}°C")
                                            
                                            # Print primer details for first few sets
                                            if len(lamp_sets) <= 3:
                                                print(f"         F3: {f3.tm:.1f}°C, B3: {b3.tm:.1f}°C")
                                                print(f"         FIP: {fip.tm:.1f}°C, BIP: {bip.tm:.1f}°C")
                                    
                                    if len(lamp_sets) >= 10000:  # Stop at 100 sets
                                        break
                                if len(lamp_sets) >= 10000:
                                    break
                            if len(lamp_sets) >= 10000:
                                break
                        if len(lamp_sets) >= 10000:
                            break
                    if len(lamp_sets) >= 10000:
                        break
                if len(lamp_sets) >= 10000:
                    break
            
            new_sets = len(lamp_sets) - initial_count
            print(f"      Found {new_sets} new sets (total: {len(lamp_sets)})")
            
            if len(lamp_sets) >= 1000:  # Stop if we have at least 10 sets
                break
        
        print(f"\n   Total tested: {combinations_tested:,} combinations")
        print(f"   Valid LAMP sets found: {len(lamp_sets)}")
        
        # Sort by Tm range (smaller is better)
        if lamp_sets:
            lamp_sets_scored = []
            for lamp_set in lamp_sets:
                tms = [lamp_set.f3.tm, lamp_set.b3.tm, lamp_set.fip.tm, lamp_set.bip.tm]
                tm_range = max(tms) - min(tms)
                score = 1.0 / (1.0 + tm_range)  # Simple scoring
                lamp_sets_scored.append((lamp_set, score, tm_range))
            
            lamp_sets_scored.sort(key=lambda x: x[1], reverse=True)
            lamp_sets = [ls for ls, _, _ in lamp_sets_scored]
            
            # Print best sets
            print(f"\n   🏆 Top 3 primer sets:")
            for i, (lamp_set, score, tm_range) in enumerate(lamp_sets_scored[:3], 1):
                print(f"   Set {i}: Tm range = {tm_range:.1f}°C, Score = {score:.3f}")
        
        return lamp_sets
    
    def emergency_primer_generation(self, sequence: str, candidates: Dict[str, List[Primer]]) -> List[LAMPPrimerSet]:
        """Emergency fallback - generate ANY valid LAMP sets regardless of balance"""
        
        print(f"\n🚨 EMERGENCY MODE: Generating primers with minimal constraints")
        
        lamp_sets = []
        
        # Just take the first valid combination we can find
        for f3 in candidates['F3'][:10]:
            for b3 in candidates['B3'][:10]:
                if f3.start >= b3.start or b3.end - f3.start < 100:
                    continue
                    
                for f1c in candidates['F1c'][:5]:
                    for f2 in candidates['F2'][:5]:
                        for b1c in candidates['B1c'][:5]:
                            for b2 in candidates['B2'][:5]:
                                
                                # Create composite primers
                                fip_seq = self.processor.reverse_complement(f1c.sequence) + f2.sequence
                                bip_seq = self.processor.reverse_complement(b1c.sequence) + b2.sequence
                                
                                try:
                                    fip_tm = self.thermo.calculate_tm(fip_seq)
                                    bip_tm = self.thermo.calculate_tm(bip_seq)
                                    fip_gc = self.thermo.calculate_gc_content(fip_seq)
                                    bip_gc = self.thermo.calculate_gc_content(bip_seq)
                                    
                                    fip = Primer(fip_seq, f1c.start, f2.end, fip_tm, fip_gc, 'FIP', 0, 0)
                                    bip = Primer(bip_seq, b1c.start, b2.end, bip_tm, bip_gc, 'BIP', 0, 0)
                                    
                                    lamp_set = LAMPPrimerSet(f3, b3, fip, bip)
                                    lamp_sets.append(lamp_set)
                                    
                                    all_tms = [f3.tm, b3.tm, fip_tm, bip_tm]
                                    tm_range = max(all_tms) - min(all_tms)
                                    
                                    print(f"   Emergency set {len(lamp_sets)}: Tm range = {tm_range:.1f}°C")
                                    print(f"      F3: {f3.sequence[:20]}... ({f3.tm:.1f}°C)")
                                    print(f"      B3: {b3.sequence[:20]}... ({b3.tm:.1f}°C)")
                                    print(f"      FIP: {fip_seq[:30]}... ({fip_tm:.1f}°C)")
                                    print(f"      BIP: {bip_seq[:30]}... ({bip_tm:.1f}°C)")
                                    
                                    if len(lamp_sets) >= 5:
                                        return lamp_sets
                                except:
                                    continue
        
        return lamp_sets
    
    def _check_very_permissive_dimers(self, primers: List[Primer]) -> bool:
        """Very permissive dimer checking - only reject EXTREMELY strong dimers"""
        very_strong_dimer_threshold = -12.0  # Only reject extremely strong dimers
        
        for i, primer1 in enumerate(primers):
            for j, primer2 in enumerate(primers[i+1:], i+1):
                dimer_dg = self.thermo.calculate_dimer_dg(primer1.sequence, primer2.sequence)
                if dimer_dg < very_strong_dimer_threshold:
                    return False  # Only reject extremely strong dimers
        return True

    def _score_and_sort_primer_sets(self, primer_sets: List[LAMPPrimerSet]) -> List[LAMPPrimerSet]:
        """Score and sort primer sets using ML or rule-based approach"""
        if self.config.use_ml_scoring and self.ml_predictor and self.ml_predictor.is_trained:
            # Use ML scoring
            scored_sets = [(lamp_set, self.ml_predictor.predict_quality(lamp_set)) 
                        for lamp_set in primer_sets]
        else:
            # Use rule-based scoring
            scored_sets = [(lamp_set, self._rule_based_score(lamp_set)) 
                        for lamp_set in primer_sets]
        
        # Sort by score (higher is better)
        scored_sets.sort(key=lambda x: x[1], reverse=True)
        
        return [lamp_set for lamp_set, score in scored_sets]

    def _rule_based_score(self, lamp_set: LAMPPrimerSet) -> float:
        """Original rule-based scoring method"""
        return self._score_primer_set(lamp_set)

    def _train_ml_model(self, primer_sets: List[LAMPPrimerSet]):
        """Train ML model on generated primer sets"""
        if not self.ml_predictor:
            return
        
        # Calculate rule-based scores for training
        scores = [self._rule_based_score(lamp_set) for lamp_set in primer_sets]
        
        # Prepare training data
        X, y = self.ml_predictor.prepare_training_data(primer_sets, scores)
        
        # Train model
        training_results = self.ml_predictor.train_model(X, y)
        
        if training_results:
            print(f"ML model trained: R² = {training_results['r2']:.3f}, "
                f"MSE = {training_results['mse']:.4f}")

    def _filter_by_specificity(self, primer_sets: List[LAMPPrimerSet]) -> Tuple[List[LAMPPrimerSet], Dict[str, List[SpecificityResult]]]:
        """Filter primer sets based on specificity"""
        specific_sets = []
        all_specificity_results = defaultdict(list)
        
        for lamp_set in primer_sets:
            set_specificity = self.specificity_checker.check_primer_set_specificity(lamp_set)
            
            # Collect results for analysis
            for primer_type, result in set_specificity.items():
                all_specificity_results[primer_type].append(result)
            
            # Check if entire set is specific
            is_set_specific = all(result.is_specific for result in set_specificity.values())
            
            if is_set_specific:
                specific_sets.append(lamp_set)
            else:
                self.filtering_stats['specificity_filtered'] += 1
        
        return specific_sets, dict(all_specificity_results)

    def _check_primer_quality(self, sequence: str) -> bool:
        """Check basic primer quality criteria"""
        # Avoid homopolymer runs > 4
        if re.search(r'(.)\1{4,}', sequence):
            return False
        
        # Avoid strong secondary structures (simplified check)
        if sequence.count('GC') > 3 or sequence.count('CG') > 3:
            return False
        
        # Check 3' end stability (no more than 2 G/C in last 5 bases)
        gc_3prime = sequence[-5:].count('G') + sequence[-5:].count('C')
        if gc_3prime > 2:
            return False
        
        return True
    
    def _check_primer_quality_flexible(self, sequence: str) -> bool:
        """More flexible primer quality criteria based on composition"""
        
        # Get sequence GC content to adjust thresholds
        gc_content = self.thermo.calculate_gc_content(sequence)
        
        # Adjust homopolymer threshold based on composition
        if gc_content < 30:  # AT-rich
            max_homopolymer = 6  # More permissive for AT-rich
        elif gc_content > 70:  # GC-rich
            max_homopolymer = 4  # Stricter for GC-rich
        else:
            max_homopolymer = 5  # Normal
        
        # Check homopolymer runs with adjusted threshold
        if re.search(rf'(.)\1{{{max_homopolymer},}}', sequence):
            return False
        
        # More flexible secondary structure check
        # Only reject if BOTH GC and CG repeats are excessive
        gc_count = sequence.count('GC')
        cg_count = sequence.count('CG')
        if gc_count > 4 and cg_count > 4:  # Both must be high
            return False
        
        # More flexible 3' end check based on composition
        if gc_content < 30:  # AT-rich sequences
            # Allow more G/C at 3' end for stability
            max_gc_3prime = 4
        else:
            max_gc_3prime = 3
        
        gc_3prime = sequence[-5:].count('G') + sequence[-5:].count('C')
        if gc_3prime > max_gc_3prime:
            return False
        
        return True

    def _check_distance_constraints(self, f1c: Primer, f2: Primer, 
                                b1c: Primer, b2: Primer) -> bool:
        """Check if primer distances meet LAMP requirements"""
        f1c_f2_dist = abs(f2.start - f1c.end)
        b1c_b2_dist = abs(b1c.start - b2.end)
        
        valid = (self.config.min_f1c_f2c_distance <= f1c_f2_dist <= self.config.max_f1c_f2c_distance and
                self.config.min_b1c_b2c_distance <= b1c_b2_dist <= self.config.max_b1c_b2c_distance)
        
        if not valid:
            self.filtering_stats['distance_filtered'] += 1
        
        return valid

    def _check_dimer_formation(self, primers: List[Primer]) -> bool:
        """Check for problematic primer-dimer formation"""
        for i, primer1 in enumerate(primers):
            for j, primer2 in enumerate(primers[i+1:], i+1):
                dimer_dg = self.thermo.calculate_dimer_dg(primer1.sequence, primer2.sequence)
                if dimer_dg < self.config.max_dimer_dg:
                    self.filtering_stats['dimer_filtered'] += 1
                    return False  # Strong dimer formation detected
        return True

    def _score_primer_set(self, lamp_set: LAMPPrimerSet) -> float:
        """Score primer set for ranking (higher is better)"""
        # Scoring based on multiple criteria
        primers = [lamp_set.f3, lamp_set.b3, lamp_set.fip, lamp_set.bip]
        
        # Tm uniformity (higher score for similar Tms)
        tms = [p.tm for p in primers]
        tm_uniformity = 1.0 / (1.0 + abs(max(tms) - min(tms)))
        
        # GC content balance
        gcs = [p.gc_content for p in primers]
        gc_balance = 1.0 / (1.0 + abs(max(gcs) - min(gcs)))
        
        # 3' end stability (prefer moderate stability, around -1.5 kcal/mol)
        dg_3primes = [p.dg_3prime for p in primers]
        stability_score = 1.0 / (1.0 + sum(abs(dg + 1.5) for dg in dg_3primes) / len(dg_3primes))
        
        # Hairpin formation penalty (less negative is better)
        dg_hairpins = [p.dg_hairpin for p in primers]
        hairpin_score = 1.0 / (1.0 + sum(abs(min(0, dg)) for dg in dg_hairpins))
        
        # Weighted composite score
        total_score = (0.3 * tm_uniformity + 
                    0.25 * gc_balance + 
                    0.25 * stability_score + 
                    0.2 * hairpin_score)
        
        return total_score
    
    def validate_primer_against_type_config(self, primer_seq: str, primer_type: LAMPPrimerType, 
                                          config: Dict) -> Tuple[bool, Dict]:
        """Validate a primer sequence against its type-specific configuration."""
        
        # Calculate primer properties
        tm = self.thermo.calculate_tm(primer_seq)
        gc = self.thermo.calculate_gc_content(primer_seq)
        length = len(primer_seq)
        dg_3prime = self.thermo.calculate_3prime_stability(primer_seq)
        dg_hairpin = self.thermo.calculate_hairpin_dg(primer_seq)
        
        properties = {
            'sequence': primer_seq,
            'tm': tm,
            'gc_content': gc,
            'length': length,
            'dg_3prime': dg_3prime,
            'dg_hairpin': dg_hairpin,
            'primer_type': primer_type.value
        }
        
        # Validate against type-specific thresholds
        validations = {
            'tm_valid': config['tm_range'][0] <= tm <= config['tm_range'][1],
            'gc_valid': config['gc_range'][0] <= gc <= config['gc_range'][1],
            'length_valid': config['length_range'][0] <= length <= config['length_range'][1],
            'dg_3prime_valid': dg_3prime >= config['dg_3prime_max'],
            'dg_hairpin_valid': dg_hairpin >= config['dg_hairpin_max']
        }
        
        is_valid = all(validations.values())
        properties['validations'] = validations
        properties['is_valid'] = is_valid
        
        # Calculate type-specific quality score
        tm_score = 1.0 - abs(tm - config['tm_optimal']) / 10.0
        gc_score = 1.0 if validations['gc_valid'] else 0.5
        length_score = 1.0 if validations['length_valid'] else 0.5
        dg_score = 1.0 if validations['dg_3prime_valid'] and validations['dg_hairpin_valid'] else 0.3
        
        quality_score = (tm_score * 0.4 + gc_score * 0.2 + length_score * 0.2 + dg_score * 0.2) * config['priority_weight']
        properties['quality_score'] = max(0.0, quality_score)
        
        return is_valid, properties
    
class TypeAwareLAMPDesigner:
    """
    Complete LAMP designer that integrates primer-type-specific thresholds
    with composition awareness and advanced filtering.
    """
    
    def __init__(self, base_config=None):
        self.base_config = base_config
        self.lamp_designer = LAMPPrimerDesigner(self.base_config)
        self.type_specific_config = PrimerTypeSpecificConfig()
        self.type_specific_designer = PrimerTypeSpecificDesigner(self.type_specific_config)
        self.thermo = None  # Initialize with your thermodynamics calculator
  
    def design_primers_with_type_awareness(self, target_sequence: str, 
                                         return_top_n: int = 5) -> List[Dict]:
        """
        Complete primer design pipeline with type-specific thresholds.
        
        Returns list of LAMP primer sets with detailed analysis.
        """
        
        print(f"🧬 TYPE-AWARE LAMP PRIMER DESIGN")
        print(f"{'='*50}")
        print(f"Target length: {len(target_sequence)} bp")
        
        # Step 1: Find candidates with type-specific thresholds
        candidates, adapted_configs = self.lamp_designer.find_primers_with_type_specific_thresholds(
            target_sequence, self.type_specific_designer, self.thermo
        )
        
        # Step 2: Create balanced LAMP sets
        lamp_sets = self.lamp_designer.create_balanced_lamp_sets_with_type_awareness(
            candidates, adapted_configs, self.thermo
        )
        
        # Step 3: Return top N sets with analysis
        top_sets = lamp_sets[:return_top_n]
        
        print(f"\n🎯 FINAL RESULTS:")
        print(f"{'='*50}")
        print(f"Generated {len(lamp_sets)} valid LAMP sets")
        print(f"Returning top {len(top_sets)} sets")
        
        # Print summary of top sets
        for i, lamp_set in enumerate(top_sets, 1):
            print(f"\nSet {i}: Tm range {lamp_set['tm_range']:.1f}°C, "
                  f"Quality {lamp_set['overall_quality']:.3f}, "
                  f"Hierarchy {lamp_set['tm_hierarchy_score']:.1f}°C")
        
        return top_sets, adapted_configs
    
    def print_lamp_set_analysis(lamp_set: Dict, adapted_configs: Dict):
        """Print detailed analysis of a LAMP primer set with type-specific validation"""
        
        print(f"\n📋 LAMP SET ANALYSIS:")
        print(f"{'='*60}")
        
        # Individual primer analysis
        primers = [
            ('F3', lamp_set['f3']),
            ('B3', lamp_set['b3']),
            ('F1c', lamp_set['f1c']),
            ('F2', lamp_set['f2']),
            ('B1c', lamp_set['b1c']),
            ('B2', lamp_set['b2'])
        ]
        
        print(f"{'Type':<6} {'Sequence':<25} {'Tm':<8} {'GC%':<6} {'Len':<4} {'Score':<6} {'Status'}")
        print(f"{'-'*6} {'-'*25} {'-'*8} {'-'*6} {'-'*4} {'-'*6} {'-'*10}")
        
        for primer_type, primer_data in primers:
            seq_display = primer_data['sequence'][:22] + "..." if len(primer_data['sequence']) > 22 else primer_data['sequence']
            status = "✅ PASS" if primer_data.get('is_valid', True) else "❌ FAIL"
            
            print(f"{primer_type:<6} {seq_display:<25} {primer_data['tm']:<8.1f} {primer_data['gc_content']:<6.1f} "
                f"{primer_data['length']:<4} {primer_data.get('quality_score', 0.0):<6.3f} {status}")
        
        # Composite primer analysis
        print(f"\n🔗 COMPOSITE PRIMERS:")
        print(f"FIP: {lamp_set['fip_seq'][:40]}... (Tm={lamp_set['fip_tm']:.1f}°C, GC={lamp_set['fip_gc']:.1f}%)")
        print(f"BIP: {lamp_set['bip_seq'][:40]}... (Tm={lamp_set['bip_tm']:.1f}°C, GC={lamp_set['bip_gc']:.1f}%)")
        
        # Tm hierarchy analysis
        print(f"\n🌡️ TM HIERARCHY ANALYSIS:")
        outer_avg = (lamp_set['f3']['tm'] + lamp_set['b3']['tm']) / 2
        inner_comp_avg = (lamp_set['f1c']['tm'] + lamp_set['b1c']['tm']) / 2
        inner_seq_avg = (lamp_set['f2']['tm'] + lamp_set['b2']['tm']) / 2
        composite_avg = (lamp_set['fip_tm'] + lamp_set['bip_tm']) / 2
        
        print(f"   Outer primers (F3/B3): {outer_avg:.1f}°C")
        print(f"   Inner complements (F1c/B1c): {inner_comp_avg:.1f}°C")
        print(f"   Inner sequences (F2/B2): {inner_seq_avg:.1f}°C")
        print(f"   Composite primers (FIP/BIP): {composite_avg:.1f}°C")
        print(f"   Hierarchy score: {lamp_set.get('tm_hierarchy_score', 0.0):.1f}°C")
    
    def analyze_primer_set(self, lamp_set: Dict, adapted_configs: Dict):
        """Provide detailed analysis of a specific primer set"""
        self.print_lamp_set_analysis(lamp_set, adapted_configs)
    
    def export_primer_set_for_synthesis(self, lamp_set: Dict) -> Dict:
        """Export primer set in format suitable for DNA synthesis ordering"""
        
        synthesis_order = {
            'F3_primer': {
                'name': 'F3_Forward_Outer',
                'sequence': lamp_set['f3']['sequence'],
                'length': lamp_set['f3']['length'],
                'tm': lamp_set['f3']['tm'],
                'gc_percent': lamp_set['f3']['gc_content']
            },
            'B3_primer': {
                'name': 'B3_Backward_Outer', 
                'sequence': lamp_set['b3']['sequence'],
                'length': lamp_set['b3']['length'],
                'tm': lamp_set['b3']['tm'],
                'gc_percent': lamp_set['b3']['gc_content']
            },
            'FIP_primer': {
                'name': 'FIP_Forward_Inner',
                'sequence': lamp_set['fip_seq'],
                'length': len(lamp_set['fip_seq']),
                'tm': lamp_set['fip_tm'],
                'gc_percent': lamp_set['fip_gc'],
                'components': f"F1c({len(lamp_set['f1c']['sequence'])}bp) + F2({len(lamp_set['f2']['sequence'])}bp)"
            },
            'BIP_primer': {
                'name': 'BIP_Backward_Inner',
                'sequence': lamp_set['bip_seq'], 
                'length': len(lamp_set['bip_seq']),
                'tm': lamp_set['bip_tm'],
                'gc_percent': lamp_set['bip_gc'],
                'components': f"B1c({len(lamp_set['b1c']['sequence'])}bp) + B2({len(lamp_set['b2']['sequence'])}bp)"
            },
            'reaction_conditions': {
                'recommended_temperature': '65°C',
                'estimated_tm_range': f"{lamp_set['tm_range']:.1f}°C",
                'primer_hierarchy': 'Outer ≥ Composite > Inner components'
            }
        }
        
        return synthesis_order
    
class ResultsWriter:
    """Handles writing primer design results to various file formats"""
    
    def __init__(self, output_dir: str):
        self.output_dir = output_dir
        os.makedirs(output_dir, exist_ok=True)
        
    def write_results(self, results: Union[List[LAMPPrimerSet], DesignResults], 
                 config: PrimerConfig, target_info: Dict,
                 full_results: Optional[List[LAMPPrimerSet]] = None):
        """Write results in multiple formats - FIXED to handle empty results"""
        
        timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
        base_name = f"lamp_primers_{timestamp}"
        
        # Determine what results we have
        if isinstance(results, DesignResults):
            primer_sets = results.primer_sets
            detailed_stats = results.get_statistics()
        else:
            primer_sets = results
            detailed_stats = None
        
        # FIXED: Handle empty results gracefully
        try:
            # Write main results (even if empty - for documentation)
            self._write_primer_sets_csv(primer_sets, f"{base_name}_top.csv", config)
            self._write_primer_sets_json(primer_sets, f"{base_name}_top.json", config, detailed_stats)
            
            # Write full results if provided
            if full_results:
                self._write_primer_sets_csv(full_results, f"{base_name}_all.csv", config)
                self._write_primer_sets_json(full_results, f"{base_name}_all.json", config)
            
            # Write detailed analysis if available
            if detailed_stats:
                self._write_analysis_report(detailed_stats, f"{base_name}_analysis.txt", config, target_info)
            
            # Write configuration
            self._write_config(config, f"{base_name}_config.json")
            
            print(f"📁 Results saved to: {self.output_dir}")
            print(f"   - Main results: {base_name}_top.csv/json")
            if full_results:
                print(f"   - All results: {base_name}_all.csv/json")
            if detailed_stats:
                print(f"   - Analysis: {base_name}_analysis.txt")
            print(f"   - Configuration: {base_name}_config.json")
            
        except Exception as e:
            print(f"⚠️  Error saving results: {e}")
            # Create minimal error report
            error_file = os.path.join(self.output_dir, f"{base_name}_error_report.txt")
            with open(error_file, 'w') as f:
                f.write(f"LAMP Primer Design Error Report\n")
                f.write(f"Generated: {datetime.now()}\n")
                f.write(f"Error: {str(e)}\n")
                f.write(f"Primer sets found: {len(primer_sets) if primer_sets else 0}\n")

    def _write_primer_sets_csv(self, primer_sets: List[LAMPPrimerSet], filename: str, config: PrimerConfig):
        """Write primer sets to CSV format"""
        filepath = os.path.join(self.output_dir, filename)
        
        with open(filepath, 'w', newline='') as f:
            writer = csv.writer(f)
            
            # Header
            headers = [
                'Set_ID', 'F3_Sequence', 'F3_Tm', 'F3_GC', 'F3_dG_3prime', 'F3_dG_hairpin',
                'B3_Sequence', 'B3_Tm', 'B3_GC', 'B3_dG_3prime', 'B3_dG_hairpin',
                'FIP_Sequence', 'FIP_Tm', 'FIP_GC', 'FIP_dG_3prime', 'FIP_dG_hairpin',
                'BIP_Sequence', 'BIP_Tm', 'BIP_GC', 'BIP_dG_3prime', 'BIP_dG_hairpin',
                'LF_Sequence', 'LF_Tm', 'LF_GC', 'LF_dG_3prime', 'LF_dG_hairpin',
                'LB_Sequence', 'LB_Tm', 'LB_GC', 'LB_dG_3prime', 'LB_dG_hairpin'
            ]
            writer.writerow(headers)
            
            # Data rows
            for i, lamp_set in enumerate(primer_sets, 1):
                row = [
                    f"Set_{i:03d}",
                    lamp_set.f3.sequence, lamp_set.f3.tm, lamp_set.f3.gc_content, 
                    lamp_set.f3.dg_3prime, lamp_set.f3.dg_hairpin,
                    lamp_set.b3.sequence, lamp_set.b3.tm, lamp_set.b3.gc_content,
                    lamp_set.b3.dg_3prime, lamp_set.b3.dg_hairpin,
                    lamp_set.fip.sequence, lamp_set.fip.tm, lamp_set.fip.gc_content,
                    lamp_set.fip.dg_3prime, lamp_set.fip.dg_hairpin,
                    lamp_set.bip.sequence, lamp_set.bip.tm, lamp_set.bip.gc_content,
                    lamp_set.bip.dg_3prime, lamp_set.bip.dg_hairpin,
                    lamp_set.lf.sequence if lamp_set.lf else "", 
                    lamp_set.lf.tm if lamp_set.lf else "",
                    lamp_set.lf.gc_content if lamp_set.lf else "",
                    lamp_set.lf.dg_3prime if lamp_set.lf else "",
                    lamp_set.lf.dg_hairpin if lamp_set.lf else "",
                    lamp_set.lb.sequence if lamp_set.lb else "",
                    lamp_set.lb.tm if lamp_set.lb else "",
                    lamp_set.lb.gc_content if lamp_set.lb else "",
                    lamp_set.lb.dg_3prime if lamp_set.lb else "",
                    lamp_set.lb.dg_hairpin if lamp_set.lb else ""
                ]
                writer.writerow(row)
    
    def _write_primer_sets_json(self, primer_sets: List[LAMPPrimerSet], filename: str, 
                           config: PrimerConfig, stats: Optional[Dict] = None):
        """Write primer sets to JSON format with proper serialization"""
        filepath = os.path.join(self.output_dir, filename)
        
        # FIXED: Convert any enum keys to strings for JSON serialization
        serializable_stats = None
        if stats:
            serializable_stats = self._make_json_serializable(stats)
        
        data = {
            "metadata": {
                "timestamp": datetime.now().isoformat(),
                "tool": "LAMP Primer Designer",
                "version": "1.0.0",
                "total_sets": len(primer_sets)
            },
            "statistics": serializable_stats,
            "primer_sets": []
        }
        
        for i, lamp_set in enumerate(primer_sets, 1):
            set_data = {
                "set_id": f"Set_{i:03d}",
                "primers": {
                    "F3": self._primer_to_dict(lamp_set.f3),
                    "B3": self._primer_to_dict(lamp_set.b3),
                    "FIP": self._primer_to_dict(lamp_set.fip),
                    "BIP": self._primer_to_dict(lamp_set.bip)
                }
            }
            
            if lamp_set.lf:
                set_data["primers"]["LF"] = self._primer_to_dict(lamp_set.lf)
            if lamp_set.lb:
                set_data["primers"]["LB"] = self._primer_to_dict(lamp_set.lb)
            
            data["primer_sets"].append(set_data)
        
        with open(filepath, 'w') as f:
            json.dump(data, f, indent=2, default=str)

    def _primer_to_dict(self, primer: Primer) -> Dict:
        """Convert primer to dictionary"""
        return {
            "sequence": primer.sequence,
            "start": primer.start,
            "end": primer.end,
            "tm": round(primer.tm, 2),
            "gc_content": round(primer.gc_content, 2),
            "dg_3prime": round(primer.dg_3prime, 3),
            "dg_hairpin": round(primer.dg_hairpin, 3),
            "primer_type": primer.primer_type
        }
    
    def _write_analysis_report(self, stats: Dict, filename: str, config: PrimerConfig, target_info: Dict):
        """Write detailed analysis report"""
        filepath = os.path.join(self.output_dir, filename)
        
        with open(filepath, 'w') as f:
            f.write("LAMP PRIMER DESIGN ANALYSIS REPORT\n")
            f.write("=" * 50 + "\n\n")
            
            f.write(f"Generated: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n")
            f.write(f"Target: {target_info.get('name', 'Unknown')}\n")
            f.write(f"Sequence Length: {target_info.get('length', 'Unknown')} bp\n\n")
            
            f.write("DESIGN PARAMETERS:\n")
            f.write("-" * 20 + "\n")
            f.write(f"Primer Length: {config.min_length}-{config.max_length} bp\n")
            f.write(f"Melting Temperature: {config.min_tm}-{config.max_tm} °C\n")
            f.write(f"GC Content: {config.min_gc}-{config.max_gc} %\n")
            f.write(f"Max ΔG (3' end): {config.max_3prime_dg} kcal/mol\n")
            f.write(f"Max ΔG (hairpin): {config.max_hairpin_dg} kcal/mol\n")
            f.write(f"Max ΔG (dimer): {config.max_dimer_dg} kcal/mol\n")
            f.write(f"ML Scoring: {'Enabled' if config.use_ml_scoring else 'Disabled'}\n")
            f.write(f"Specificity Check: {'Enabled' if config.enable_specificity_check else 'Disabled'}\n\n")
            
            f.write("DESIGN STATISTICS:\n")
            f.write("-" * 20 + "\n")
            f.write(f"Total Candidates Evaluated: {stats['total_candidates']}\n")
            f.write(f"Total Combinations Tested: {stats['total_combinations']}\n")
            f.write(f"Valid Primer Sets Found: {stats['valid_sets']}\n")
            f.write(f"Success Rate: {(stats['valid_sets']/stats['total_combinations']*100):.2f}%\n\n")
            
            f.write("FILTERING BREAKDOWN:\n")
            f.write("-" * 20 + "\n")
            filtering = stats['filtering_stats']
            for filter_type, count in filtering.items():
                f.write(f"{filter_type.replace('_', ' ').title()}: {count}\n")
            
            if stats['score_range'][0] != stats['score_range'][1]:
                f.write(f"\nScore Range: {stats['score_range'][0]:.3f} - {stats['score_range'][1]:.3f}\n")
            
            if 'specificity_stats' in stats:
                f.write("\nSPECIFICITY ANALYSIS:\n")
                f.write("-" * 20 + "\n")
                for primer_type, spec_stats in stats['specificity_stats'].items():
                    f.write(f"{primer_type}: {spec_stats['specific_primers']}/{spec_stats['total_checked']} specific\n")
    
    def _write_config(self, config: PrimerConfig, filename: str):
        """Write configuration to JSON file"""
        filepath = os.path.join(self.output_dir, filename)
        
        config_dict = {
            "primer_constraints": {
                "min_length": config.min_length,
                "max_length": config.max_length,
                "min_tm": config.min_tm,
                "max_tm": config.max_tm,
                "min_gc": config.min_gc,
                "max_gc": config.max_gc
            },
            "thermodynamic_constraints": {
                "max_3prime_dg": config.max_3prime_dg,
                "max_hairpin_dg": config.max_hairpin_dg,
                "max_dimer_dg": config.max_dimer_dg
            },
            "reaction_conditions": {
                "na_conc": config.na_conc,
                "k_conc": config.k_conc,
                "mg_conc": config.mg_conc,
                "temperature": config.temperature
            },
            "advanced_options": {
                "use_ml_scoring": config.use_ml_scoring,
                "enable_specificity_check": config.enable_specificity_check,
                "sequence_type": config.seq_type.value
            }
        }
        
        with open(filepath, 'w') as f:
            json.dump(config_dict, f, indent=2)

    def _make_json_serializable(self, obj):
        """Convert objects to JSON-serializable format"""
        if isinstance(obj, dict):
            new_dict = {}
            for key, value in obj.items():
                # Convert LAMPPrimerType and other enums to strings
                if hasattr(key, 'value'):  # It's an enum
                    str_key = str(key.value)
                else:
                    str_key = str(key)
                new_dict[str_key] = self._make_json_serializable(value)
            return new_dict
        elif isinstance(obj, list):
            return [self._make_json_serializable(item) for item in obj]
        elif hasattr(obj, 'value'):  # It's an enum
            return obj.value
        else:
            return obj

def parse_args():
    """Parse command line arguments"""
    parser = ArgumentParser(description='LAMP Primer Designer - Complete Tool')
    
    # Input/Output
    parser.add_argument('-i', '--input', required=True, 
                       help='Input sequence (FASTA file or string)')
    parser.add_argument('-o', '--output', required=True, 
                       help='Output directory for results')
    
    # Output options
    parser.add_argument('-f', '--format', default='top_n', 
                       choices=['top_n', 'detailed'], 
                       help='Result format: top_n or detailed (default: top_n)')
    parser.add_argument('--full', action='store_true', 
                       help='Save ALL primer combinations (disables quality filtering)')
    parser.add_argument('-n', '--top_n', type=int, default=3,
                       help='Number of top primer sets to return (default: 3)')
    
    # UPDATED: More realistic primer constraints based on actual LAMP primer data
    parser.add_argument('-l', '--min_length', type=int, default=15, 
                       help='Minimum primer length (default: 15)')
    parser.add_argument('-L', '--max_length', type=int, default=40, 
                       help='Maximum primer length (default: 40)')
    parser.add_argument('-t', '--min_tm', type=float, default=50.0, 
                       help='Minimum melting temperature (default: 50.0)')
    parser.add_argument('-T', '--max_tm', type=float, default=80.0, 
                       help='Maximum melting temperature (default: 80.0)')
    
    # FIXED: Realistic GC content range with both min and max
    parser.add_argument('-c', '--min_gc', type=float, default=20.0, 
                       help='Minimum primer GC %% (default: 20.0)')
    parser.add_argument('-C', '--max_gc', type=float, default=80.0, 
                       help='Maximum primer GC %% (default: 80.0)')
    
    # Thermodynamic constraints (more permissive)
    parser.add_argument('-g', '--max_3prime_dg', type=float, default=-0.5, 
                       help='Maximum ΔG for 3\' end stability (default: -0.5)')
    parser.add_argument('-G', '--max_hairpin_dg', type=float, default=-1.0, 
                       help='Maximum ΔG for hairpin formation (default: -1.0)')
    parser.add_argument('-d', '--max_dimer_dg', type=float, default=-3.0, 
                       help='Maximum ΔG for primer-dimer formation (default: -3.0)')
    
    # Advanced options
    parser.add_argument('--ml_scoring', action='store_true', 
                       help='Enable machine learning scoring')
    parser.add_argument('-s', '--seq_type', default='normal', 
                       choices=['normal', 'bisulfite', 'pyrosequencing'], 
                       help='Sequencing type (default: normal)')
    parser.add_argument('--loop_primers', action='store_true', 
                       help='Design loop primers (LF/LB)')
    parser.add_argument('--no_specificity', action='store_true', 
                       help='Disable specificity checking')
    
    # Type-specific primer design
    
    parser.add_argument('--enable-type-specific', action='store_true',
                       help='Enable primer-type-specific threshold optimization')
    
    parser.add_argument('--tm-hierarchy-weight', type=float, default=0.1,
                       help='Weight for Tm hierarchy scoring (default: 0.1)')
    
    parser.add_argument('--outer-tm-boost', type=float, default=0.0,
                       help='Additional Tm boost for outer primers (default: 0.0)')
    
    parser.add_argument('--composite-length-max', type=int, default=50,
                       help='Maximum length for composite primers (default: 50)')
    
    parser.add_argument('--print-type-analysis', action='store_true',
                       help='Print detailed type-specific analysis for each primer')
    
    return parser.parse_args()

def parse_fasta(filepath: str) -> Tuple[str, str]:
    """Parse FASTA file and return header and sequence"""
    header = ""
    sequence = ""
    
    with open(filepath, 'r') as f:
        for line in f:
            line = line.strip()
            if line.startswith('>'):
                header = line[1:]  # Remove '>' character
            elif line:
                sequence += line.upper()
    
    return header, sequence

def integrate_type_specific_design(existing_lamp_designer, target_sequence: str):
    """Integration function to add primer-type-specific design to existing LAMPPrimerDesigner."""
    
    print("🎯 Enhancing LAMP designer with primer-type-specific thresholds...")
    
    # Create type-specific configuration
    type_specific_config = PrimerTypeSpecificConfig()
    type_specific_designer = PrimerTypeSpecificDesigner(type_specific_config)
    
    # Find primers with type-specific thresholds
    candidates, adapted_configs = existing_lamp_designer.find_primers_with_type_specific_thresholds(
        target_sequence, type_specific_designer
    )
    
    # Create balanced LAMP sets
    lamp_sets = existing_lamp_designer.create_balanced_lamp_sets_with_type_awareness(
        candidates, adapted_configs
    )
    
    return lamp_sets, adapted_configs

def main():
    """Main function for LAMP primer design"""
    print("🧬 === LAMP Primer Designer - Complete Tool ===\n")
    
    args = parse_args()
    
    # Parse input
    if os.path.isfile(args.input):
        header, target = parse_fasta(args.input)
        target_info = {"name": header or os.path.basename(args.input), "length": len(target)}
        print(f"📁 Loaded sequence from file: {args.input}")
        print(f"   - Name: {target_info['name']}")
        print(f"   - Length: {target_info['length']} bp")
    else:
        target = args.input.upper()
        target_info = {"name": "Direct input", "length": len(target)}
        print(f"📝 Using direct sequence input ({len(target)} bp)")
    
    # Parse sequence type
    seq_type_map = {
        'normal': SequenceType.NORMAL,
        'bisulfite': SequenceType.BISULFITE,
        'pyrosequencing': SequenceType.PYROSEQUENCING
    }
    seq_type = seq_type_map[args.seq_type]
    
    # Create configuration
    config = PrimerConfig(
        min_length=args.min_length,
        max_length=args.max_length,
        min_tm=args.min_tm,
        max_tm=args.max_tm,
        min_gc=float(args.min_gc),      
        max_gc=float(args.max_gc),      
        max_3prime_dg=args.max_3prime_dg,
        max_hairpin_dg=args.max_hairpin_dg,
        max_dimer_dg=args.max_dimer_dg,
        use_ml_scoring=args.ml_scoring,
        enable_specificity_check=not args.no_specificity,
        seq_type=seq_type
    )
    
    print(f"⚙️  Configuration:")
    print(f"   - Primer length: {config.min_length}-{config.max_length} bp")
    print(f"   - Melting temp: {config.min_tm}-{config.max_tm} °C")
    print(f"   - ML scoring: {'Enabled' if config.use_ml_scoring else 'Disabled'}")
    print(f"   - Specificity check: {'Enabled' if config.enable_specificity_check else 'Disabled'}")
    print(f"   - Sequence type: {seq_type.value}")
    print(f"   - Quality filtering: {'Disabled (--full mode)' if args.full else 'Enabled'}")
    
    # Create designer
    print(f"\n🔬 Starting primer design...")
    designer = LAMPPrimerDesigner(config)
    
    # Modify designer behavior for --full mode
    if args.full:
        print(f"📊 Full mode: Generating ALL primer combinations (no filtering)...")
        designer.disable_filtering = True
    
    # === MAIN DESIGN LOGIC ===
    # Choose design method based on arguments
    if hasattr(args, 'enable_type_specific') and args.enable_type_specific:
        print("🎯 Using primer-type-specific thresholds...")
        
        try:
            # Enhanced type-aware design with existing designer
            lamp_sets, adapted_configs = integrate_type_specific_design(designer, target)
            
            # Convert to standard format - FIXED VERSION
            if args.format == 'detailed':
                # Proper calculation of total_candidates
                total_candidates = 0
                if hasattr(designer, 'filtering_stats'):
                    total_candidates = designer.filtering_stats.get('total_candidates', 0)
                elif lamp_sets:
                    # Estimate from lamp_sets if filtering_stats not available
                    total_candidates = len(lamp_sets) * 10  # rough estimate
                
                # Convert lamp_sets to LAMPPrimerSet objects first
                converted_lamp_sets = convert_type_aware_to_lamp_sets(lamp_sets[:args.top_n]) if lamp_sets else []
                
                # Create DesignResults object for detailed format
                results = DesignResults(
                    primer_sets=converted_lamp_sets,
                    total_candidates=total_candidates,
                    total_combinations=len(lamp_sets) if lamp_sets else 0,
                    filtering_stats=designer.filtering_stats.copy() if hasattr(designer, 'filtering_stats') else {},
                    scores=[lamp_set.get('overall_quality', 0.0) for lamp_set in lamp_sets[:args.top_n]] if lamp_sets else [],
                    composition_analysis=adapted_configs,  # Add directly in constructor
                    design_method="type_specific"  # Add directly in constructor
                )
                
            else:
                # Convert dict format to LAMPPrimerSet objects
                results = convert_type_aware_to_lamp_sets(lamp_sets[:args.top_n]) if lamp_sets else []
            
            # Store full results for --full mode
            full_results = convert_type_aware_to_lamp_sets(lamp_sets) if args.full and lamp_sets else None
            
            # Print detailed analysis if requested
            if hasattr(args, 'print_type_analysis') and args.print_type_analysis and lamp_sets:
                print_lamp_set_analysis(lamp_sets[0], adapted_configs)
            
            # Export for synthesis
            if lamp_sets:
                os.makedirs(args.output, exist_ok=True)
                synthesis_order = export_lamp_set_for_synthesis(lamp_sets[0])
                synthesis_file = os.path.join(args.output, "synthesis_order.json")
                with open(synthesis_file, 'w') as f:
                    json.dump(synthesis_order, f, indent=2)
                print(f"💾 Synthesis order saved to: {synthesis_file}")
            
            print(f"✅ Type-specific design complete - {len(lamp_sets) if lamp_sets else 0} sets found")
            
        except Exception as e:
            print(f"❌ Error in type-specific design: {e}")
            print("🔄 Falling back to standard design method...")
            # Fall back to standard method
            results = designer.design_primers(target, ResultFormat.TOP_N, top_n=args.top_n)
            full_results = None
    
    else:
        # === STANDARD DESIGN METHOD ===
        print("🔬 Using standard primer design...")
        
        try:
            # Get results based on format
            if args.format == 'detailed':
                results = designer.design_primers(target, ResultFormat.DETAILED)
                print(f"✅ Detailed analysis complete")
                
                # Get all results for --full mode
                full_results = results.primer_sets if args.full else None
                
            else:
                results = designer.design_primers(target, ResultFormat.TOP_N, top_n=args.top_n)
                print(f"✅ Found top {len(results) if isinstance(results, list) else 0} primer sets")
                
                # Get all results for --full mode
                if args.full:
                    full_results = designer.design_primers(target, ResultFormat.ALL)
                    if isinstance(full_results, DesignResults):
                        full_results = full_results.primer_sets
                    print(f"   - Total primer combinations: {len(full_results) if full_results else 0}")
                else:
                    full_results = None
                    
        except Exception as e:
            print(f"❌ Error in standard design: {e}")
            # Return empty results
            if args.format == 'detailed':
                results = DesignResults(
                    primer_sets=[],
                    total_candidates=0,
                    total_combinations=0,
                    filtering_stats={},
                    scores=[]
                )
            else:
                results = []
            full_results = None

    # === ADD LOOP PRIMERS (if requested) ===
    if hasattr(args, 'loop_primers') and args.loop_primers:
        print(f"🔄 Adding loop primers...")
        loop_designer = LoopPrimerDesigner(config)
        
        # Add loop primers to main results
        if isinstance(results, DesignResults):
            enhanced_sets = []
            for lamp_set in results.primer_sets:
                lf, lb = loop_designer.design_loop_primers(target, lamp_set)
                enhanced_set = LAMPPrimerSet(
                    lamp_set.f3, lamp_set.b3, lamp_set.fip, lamp_set.bip, lf, lb
                )
                enhanced_sets.append(enhanced_set)
            results.primer_sets = enhanced_sets
            
        elif isinstance(results, list):
            enhanced_sets = []
            for lamp_set in results:
                lf, lb = loop_designer.design_loop_primers(target, lamp_set)
                enhanced_set = LAMPPrimerSet(
                    lamp_set.f3, lamp_set.b3, lamp_set.fip, lamp_set.bip, lf, lb
                )
                enhanced_sets.append(enhanced_set)
            results = enhanced_sets
        
        # Add loop primers to full results
        if full_results:
            enhanced_full = []
            for lamp_set in full_results:
                lf, lb = loop_designer.design_loop_primers(target, lamp_set)
                enhanced_set = LAMPPrimerSet(
                    lamp_set.f3, lamp_set.b3, lamp_set.fip, lamp_set.bip, lf, lb
                )
                enhanced_full.append(enhanced_set)
            full_results = enhanced_full
        
        print(f"   - Loop primers added to all sets")
    
    # === SAVE RESULTS ===
    print(f"\n💾 Saving results...")
    
    # Ensure output directory exists
    os.makedirs(args.output, exist_ok=True)
    
    writer = ResultsWriter(args.output)
    writer.write_results(results, config, target_info, full_results)
    
    # === DISPLAY SUMMARY ===
    print(f"\n🎯 Summary:")
    if isinstance(results, DesignResults):
        stats = results.get_statistics()
        print(f"   - Valid primer sets found: {stats['valid_sets']}")
        print(f"   - Candidates evaluated: {stats['total_candidates']}")
        if stats['total_combinations'] > 0:
            success_rate = (stats['valid_sets']/stats['total_combinations']*100)
            print(f"   - Success rate: {success_rate:.2f}%")
        
        # Show type-specific info if available
        if hasattr(results, 'design_method') and results.design_method == "type_specific":
            print(f"   - Design method: Type-specific thresholds ✨")
            
    else:
        print(f"   - Top primer sets: {len(results) if results else 0}")
    
    if full_results:
        print(f"   - Total valid sets: {len(full_results)}")
    
    print(f"\n🧬 LAMP primer design complete!")
    print(f"📁 Results saved to: {args.output}")
    
    # === RETURN RESULTS FOR TESTING ===
    return results, full_results


def convert_type_aware_to_lamp_sets(lamp_dicts: List[Dict]) -> List[LAMPPrimerSet]:
    """Convert type-aware dict format to LAMPPrimerSet objects"""
    lamp_sets = []
    
    for lamp_dict in lamp_dicts:
        # Create Primer objects from dict data
        f3 = Primer(
            sequence=lamp_dict['f3']['sequence'],
            start=0,  # Position data may not be available in dict format
            end=len(lamp_dict['f3']['sequence']),
            tm=lamp_dict['f3']['tm'],
            gc_content=lamp_dict['f3']['gc_content'],
            primer_type='F3',
            dg_3prime=lamp_dict['f3'].get('dg_3prime', 0.0),
            dg_hairpin=lamp_dict['f3'].get('dg_hairpin', 0.0)
        )
        
        b3 = Primer(
            sequence=lamp_dict['b3']['sequence'],
            start=0,
            end=len(lamp_dict['b3']['sequence']),
            tm=lamp_dict['b3']['tm'],
            gc_content=lamp_dict['b3']['gc_content'],
            primer_type='B3',
            dg_3prime=lamp_dict['b3'].get('dg_3prime', 0.0),
            dg_hairpin=lamp_dict['b3'].get('dg_hairpin', 0.0)
        )
        
        fip = Primer(
            sequence=lamp_dict['fip_seq'],
            start=0,
            end=len(lamp_dict['fip_seq']),
            tm=lamp_dict['fip_tm'],
            gc_content=lamp_dict['fip_gc'],
            primer_type='FIP',
            dg_3prime=0.0,  # Could calculate if needed
            dg_hairpin=0.0
        )
        
        bip = Primer(
            sequence=lamp_dict['bip_seq'],
            start=0,
            end=len(lamp_dict['bip_seq']),
            tm=lamp_dict['bip_tm'],
            gc_content=lamp_dict['bip_gc'],
            primer_type='BIP',
            dg_3prime=0.0,
            dg_hairpin=0.0
        )
        
        lamp_sets.append(LAMPPrimerSet(f3, b3, fip, bip))
    
    return lamp_sets


def export_lamp_set_for_synthesis(lamp_dict: Dict) -> Dict:
    """Export type-aware LAMP set for DNA synthesis ordering"""
    
    synthesis_order = {
        'F3_primer': {
            'name': 'F3_Forward_Outer',
            'sequence': lamp_dict['f3']['sequence'],
            'length': lamp_dict['f3']['length'],
            'tm': round(lamp_dict['f3']['tm'], 1),
            'gc_percent': round(lamp_dict['f3']['gc_content'], 1),
            'primer_type': 'Outer'
        },
        'B3_primer': {
            'name': 'B3_Backward_Outer', 
            'sequence': lamp_dict['b3']['sequence'],
            'length': lamp_dict['b3']['length'],
            'tm': round(lamp_dict['b3']['tm'], 1),
            'gc_percent': round(lamp_dict['b3']['gc_content'], 1),
            'primer_type': 'Outer'
        },
        'FIP_primer': {
            'name': 'FIP_Forward_Inner',
            'sequence': lamp_dict['fip_seq'],
            'length': len(lamp_dict['fip_seq']),
            'tm': round(lamp_dict['fip_tm'], 1),
            'gc_percent': round(lamp_dict['fip_gc'], 1),
            'primer_type': 'Inner_Composite',
            'components': f"F1c({len(lamp_dict['f1c']['sequence'])}bp) + F2({len(lamp_dict['f2']['sequence'])}bp)"
        },
        'BIP_primer': {
            'name': 'BIP_Backward_Inner',
            'sequence': lamp_dict['bip_seq'], 
            'length': len(lamp_dict['bip_seq']),
            'tm': round(lamp_dict['bip_tm'], 1),
            'gc_percent': round(lamp_dict['bip_gc'], 1),
            'primer_type': 'Inner_Composite',
            'components': f"B1c({len(lamp_dict['b1c']['sequence'])}bp) + B2({len(lamp_dict['b2']['sequence'])}bp)"
        },
        'reaction_conditions': {
            'recommended_temperature': '65°C',
            'estimated_tm_range': f"{lamp_dict['tm_range']:.1f}°C",
            'primer_hierarchy': 'Outer ≥ Composite > Inner components',
            'design_method': 'Type-specific thresholds'
        },
        'quality_metrics': {
            'tm_hierarchy_score': round(lamp_dict.get('tm_hierarchy_score', 0.0), 2),
            'overall_quality': round(lamp_dict.get('overall_quality', 0.0), 3),
            'tm_balance': f"±{lamp_dict['tm_range']/2:.1f}°C"
        }
    }
    
    return synthesis_order

if __name__ == "__main__":
    main()
