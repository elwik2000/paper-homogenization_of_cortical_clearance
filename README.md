# Supplementary code: Well-posedness and homogenization of a waste-protein clearance model in cortical tissue

This repository contains supplementary code for the paper
> Bulíček, M., Jüngel, A., Rognes, M. E. & Wintrová, L.
> Well-posedness and homogenization of a waste-protein clearance model in cortical tissue.


## Abstract
We formulate and analyze a mathematical model for the clearance of waste proteins from cortical tissue. The model tracks three species of waste proteins, distinguished by their size, and incorporates constant production at neuronal surfaces, reaction–diffusion in the bulk of interstitial fluid, and nonlinear dynamics on capillary surfaces. We establish the existence, boundedness, and uniqueness of a solution to the microscopic model. 
Using two-scale convergence, we derive a homogenized model and prove that the microscopic solutions converge to the corresponding macroscopic limit. The resulting model reveals how the microscopic arrangement and relative volume fractions of neurons and capillaries affect waste-protein transport and clearance at the tissue scale. Numerical simulations illustrate these effects and demonstrate how variations in cortical microstructure translate into differences in the macroscopic clearance dynamics.

## Getting started

First install `uv`
```bash
curl -LsSf https://astral.sh/uv/install.sh | sh
```
Then create a virtual environment and install dependencies

```bash
uv venv
source .venv/bin/activate
uv sync --locked --all-extras --dev
```

We also provide a pre-build Docker image which can be used to run the the code in this repository. First thing you need to do is to ensure that you have [docker installed](https://docs.docker.com/get-docker/).

To start an interactive docker container you can execute the following command

```bash
docker run --rm -it ghcr.io/scientificcomputing/example-paper:latest
```

### Pre-processing
Add steps for pre-processing, e.g

```
cd code
python3 pre-processing.py
```

### Running simulation
Add steps for running simulations, e.g

```
cd code
python3 run_all.py
```


### Postprocessing
Add steps for postprocessing / reproducing figures and tables in the paper, e.g

```
cd code
python3 postprocess.py
```

## Citation

```
@software{Lisa_My_Research_Software_2017,
  author = {Lisa, Mona and Bot, Hew},
  doi = {10.5281/zenodo.1234},
  month = {12},
  title = {{My Research Software}},
  url = {https://github.com/scientificcomputing/example-paper},
  version = {2.0.4},
  year = {2017}
}
```


## Having issues
If you have any troubles please file and issue in the GitHub repository.

## License
MIT
