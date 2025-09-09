## Luo-Rudy Finitewave model

The Luo-Rudy-91 model is a mathematical representation of cardiac action potentials, designed to simulate the electrical activity of heart cells. This model captures the dynamics of ion currents and their interactions, providing insights into the mechanisms of cardiac excitability and conduction. It is particularly useful for studying arrhythmias and the effects of various pharmacological agents on cardiac function.

This model implementation can be used separately from the Finitewave, allowing for standalone simulations and testing of the Luo-Rudy-91 dynamics without the need for the entire framework.

### Reference
Luo, C. H., & Rudy, Y. (1991). A model of the ventricular cardiac action potential.
Depolarization, repolarization, and their interaction. Circulation Research, 68(6), 1501-1526.

DOI: https://doi.org/10.1161/01.res.68.6.1501

### How to use
```bash
python -m examples.luo_rudy_91_example
```

This will initiate a simulation using the Luo-Rudy91 model and display the results.

### How to test
```bash
python -m pytest -q
```

### Repository structure
```text
.
├── luo_rudy_91/                  # Luo-Rudy-91 model equations package (ops.py)
│   ├── __init__.py
│   └── ops.py                       
├── implementation/               # 0D model implementation
│   ├── __init__.py
│   └── luo_rudy_91_0d.py
├── example/
│   └── luo_rudy_91_example.py    # minimal script to run a short trace
├── tests/
│   └── luo_rudy_91_test.py       # Luo-Rudy-91 model test
├── .gitignore
├── LICENSE                       # MIT
├── pyproject.toml                   
└── README.md                     # this file
```

### Variables
Model state variables: description, units and ranges (optional)
- `u`   - Transmembrane potential (mV)
- `m`   - Activation variable for sodium channels
- `h`   - Inactivation variable for sodium channels
- `j`   - Inactivation variable for sodium channels
- `d`   - Activation variable for calcium channels
- `f`   - Inactivation variable for calcium channels
- `x`   - Activation variable for potassium channels
- `cai` - Intracellular calcium concentration (mM)

### Parameters
Ion Channel Conductances (mS/µF)
- `gna = 23.0`    - Fast sodium (Na+) conductance
- `gsi = 0.09`    - Slow inward calcium (Ca2+) conductance
- `gk  = 0.282`   - Time-dependent potassium (K+) conductance
- `gk1 = 0.6047`  - Inward rectifier potassium (K1) conductance
- `gkp = 0.0183`  - Plateau potassium (Kp) conductance
- `gb  = 0.03921` - Background conductance (leak current)

Extracellular and Intracellular Ion Concentrations (mM)
- `ko  = 5.4`   - Extracellular potassium concentration
- `ki  = 145.0` - Intracellular potassium concentration
- `nai = 18.0`  - Intracellular sodium concentration
- `nao = 140.0` - Extracellular sodium concentration
- `cao = 1.8`   - Extracellular calcium concentration

Physical Constants
- `R = 8.314` - Universal gas constant (J/(mol·K))
- `T = 310.0` - Temperature (Kelvin, 37°C)
- `F = 96.5`  - Faraday constant (C/mmol)

Ion Permeability Ratios
- `PR_NaK = 0.01833` - Na+/K+ permeability ratio

Equilibrium potentials (mV)
- `E_Na = (R * T / F) * log(nao / nai)`
- `E_K1 = (R * T / F) * log(ko / ki)` 

