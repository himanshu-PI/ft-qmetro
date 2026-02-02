Python scripts for [Achieving the Heisenberg limit using fault-tolerant quantum error correction](https://arxiv.org/abs/2601.05457)

## Installation
This repository uses [Conda](https://www.anaconda.com/docs/getting-started/miniconda/main) to manage dependencies.

```bash
git clone https://github.com/himanshu-PI/ft-qmetro.git
cd ft-qmetro
conda env create -f env.yml
conda activate ft-qmetro
```

## Repository Structure
```
.
├── fi/                     # Fisher Information Calculation
│   ├── data/
│   │   └── calc/           # Presaved data 
│   │       └── .txt
│   ├── essential.py        # helper functions
│   ├── fi.py               # calculates FI
│   └── mean-sq-phase.py    # calculates average mean square magnetization for FI
├── state-prep/
│   ├── correlated-noise/   # case study : correlated noise
│   │   ├── data/       
│   │   │   └── .txt
│   │   ├── funcs.py        # helper functions
│   │   └── main.py         # calcualtes average mean square magnetization
│   ├── data/
│   │   └── calc/`          # presaved data
│   │       └── .txt
│   ├── essential.py        # helper functions
│   └── mean-sq-phase.py    # calcualtes average mean square magnetization
├── CITATION.cff
├── LICENSE
├── env.yml
└── README.md
```


## Usage

- Calculate average mean square magnetization 
    ```bash
    cd state-prep
    python mean-sq-phase.py
    ```
- Calculate Fisher information
    ```bash
    cd fi
    python mean-sq-phase.py
    python fi.py
    ```

