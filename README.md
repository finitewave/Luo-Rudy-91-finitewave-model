## Finitewave model template (replace with the model name)

Add model description here and fill the sections below.

This model implementation can be used separately from the Finitewave, allowing for standalone simulations and testing of the model dynamics without the need for the entire framework.


### Reference
Paper, Authors, DOI.

### How to use (quickstart)
```bash
python -m example/model_example
```

### How to test
```bash
python -m pytest -q
```

### Repository structure
```text
.
├── model_template/                  # equations package (ops.py)
│   ├── __init__.py
│   └── ops.py                       # fill with the model equations (pure functions)
├── implementation/                  # 0D model implementation
│   ├── __init__.py
│   └── model_0d.py
├── example/
│   └── model_example.py             # minimal script to run a short trace
├── tests/
│   └── test.py                      # smoke test; extend with reproducibility checks
├── .gitignore
├── LICENSE                          # MIT
├── pyproject.toml                   # placeholders to replace
└── README.md                        # this file
```

### Variables
Model state variables: description, units and ranges (optional)
- `u` — ...

### Parameters
Parameters and their defualt values
- `par` - ...

### Model Contributor TODO (template repository only)

In pyproject.toml, replace "model_template" with the actual model id. 
It must match the name of the directory where ops.py is located.
```toml
[project.entry-points."finitewave.models"]
model_template = "finitewave_models.model_template"
```
