# NSD2 Paper at Nature

A reproducible pipeline to load, process, and visualize NSD2 gene expression and protein activity data for the Nature paper submission.

![Diffusion Map Leiden](assets/diffmap_scRNASeq.png)

## Repository Structure

```
nsd2-paper-at-nature/
├── src/
│   ├── load_dataset.py          # Script to load raw data, process with custom functions, and save outputs
│   ├── processing_functions.py  # Helper functions for data processing and formatting
│   └── print_charts.ipynb       # Jupyter notebook to generate figures for the manuscript
└── README.md                    # This file
```

## Prerequisites

- Python 3.8+
- Required Python packages listed at the top of `src/load_dataset.py` and `src/processing_functions.py` (e.g., `scanpy`, `pandas`, `numpy`, etc.)
- A Jupyter environment or [JupyterLab](https://jupyter.org/) to run the notebook

## Usage

Follow these two main steps to reproduce the data files and figures:

### 1. Load and Process Data

1. Open `src/load_dataset.py` in a text editor.
2. Locate the `new_data_dir` variable and set it to your data folder, e.g.:
   ```python
   new_data_dir = "/path/to/your/data/folder"
   ```
3. Run the script to download/process raw data and save results:
   ```bash
   python src/load_dataset.py
   ```
4. After completion, you will find new subdirectories under `new_data_dir`. These will be used to generate the figures.

### 2. Generate Figures

1. Launch Jupyter Notebook:
   ```bash
   jupyter notebook src/print_charts.ipynb
   ```
2. In the first cell of `print_charts.ipynb`, update the `base_path` variable to the same `new_data_dir` used above:
   ```python
   base_path = "/path/to/your/data/folder"
   ```
3. Execute all cells to recreate the figures for the paper. Outputs (plots) will be saved and rendered inline.

## File Descriptions

- **load\_dataset.py**:

  - Downloads or reads raw gene expression datasets.
  - Uses functions from `processing_functions.py` to clean, filter, and annotate data.

- **processing\_functions.py**:

  - Contains functions for:
    - Reading raw files
    - Generating and saving Gene Expression Data
    - Generating and saving Protein Activity Data

- **print\_charts.ipynb**:

  - Reads processed data from `base_path`.
  - Generates manuscript figures (e.g., DiffMaps, heatmaps, Dotplots).
