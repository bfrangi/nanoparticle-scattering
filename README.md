# Nanoparticle Scattering Simulator

Code for simulating Mie scattering by nanoparticles.

## Quick Start

To try this project out, follow these steps:

1. Clone the repo
    ```
    git clone https://github.com/bfrangi/nanoparticle-scattering.git
    ```
2. Navigate to the project folder: 
    ```
    cd nanoparticle-scattering
    ```
3. Install `pip`. On most Linux systems, this can be done with:
    ```
    sudo apt install python3-pip
    ```
4. Install `virtualenv` and create a virtual environment:
    ```
    sudo apt install python3-venv
    python3 -m venv .venv
    source .venv/bin/activate
    ```
5. Install dependencies:
    ```
    pip install -r requirements.txt
    ```
6. [Optional] If you want to use LaTeX for plotting, install the `texlive` package:

   ```bash
   sudo apt-get update
   sudo apt-get install texlive-latex-extra texlive-fonts-extra dvipng
   ```

   **Note:** Using LaTeX for plotting is optional, and slows down the plotting process significantly.
   It is recommended only for generating publication-quality figures.


## Usage

The main simulation scripts are located in the `src/` folder. You can run them using Python. For example:
```
python src/si-intensity-interference.py
```

These scripts will generate plots and save them in the `figures/` folder.