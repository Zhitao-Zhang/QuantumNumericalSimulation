# CIMSim - Coherent Ising Machine Simulator

## Project Structure

```
CIMSim/
├── cimsim/                  # Our CIM simulator implementation (TBD)
├── baseline/                # SOTA baseline: mcmahon-lab/cim-optimizer
│   └── cim_optimizer/
│       ├── AHC.py           # Amplitude-Heterogeneity Correction
│       ├── CAC.py           # Chaotic Amplitude Control
│       ├── extAHC.py        # AHC with external field
│       └── CIM_helper.py    # Utility functions
├── benchmark/
│   ├── instances.py         # Instance generators (random, MAX-CUT, SK model)
│   └── evaluate.py          # Evaluation framework
├── instances/               # Generated benchmark instances (.npz)
└── requirements.txt
```

## Quick Start

```bash
# Generate benchmark instances
python3 benchmark/instances.py

# Run baseline evaluation
python3 benchmark/evaluate.py
```

## Baseline: mcmahon-lab/cim-optimizer

The baseline implements three CIM algorithm variants from McMahon Lab (Cornell):

- **CAC** (Chaotic Amplitude Control) - Leleu et al., Commun Phys 2021
- **AHC** (Amplitude-Heterogeneity Correction) - Leleu et al., PRL 2019
- **extAHC** - AHC with external field support

## Benchmark Instances

15 instances across three categories:
- **Random Ising**: N = 20, 50, 100, 200, 500, 1000
- **MAX-CUT**: N = 20, 50, 100, 200, 500
- **SK Model**: N = 50, 100, 200, 500
