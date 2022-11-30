# GIRFReco.jl: An open-source pipeline for spiral MRI Reconstruction in Julia

### 0.7.1 Requirements
* `dev MRIReco`
* `dev MRIFiles`
* `dev MRIBase`
* `add ImageUtils`
* `add MRICoilSensitivities`
* `add FileIO`
* ... more to come when debugging works 

This package provides an image reconstruction pipeline for real-world MRI use cases, such as spiral diffusion imaging. It is completely implemented in Julia using original code and external packages, e.g., [MRIReco.jl](https://magneticresonanceimaging.github.io/MRIReco.jl/latest/) for the main iterative reconstruction tasks and [MRIGradients.jl](https://github.com/BRAIN-TO/MRIGradients.jl) for the corrections of system imperfections via the Gradient Impulse Response Function (GIRF).

This repository includes a working example for spiral reconstruction with GIRF correction of trajectory (kx-ky-kz, dubbed as k1) and B0 eddy currents (k0), iterative reconstruction (cg-SENSE) and a Cartesian reconstruction example for coil sensitivity and off-resonance map calculation.

The data for the phantom reconstruction (`SpiralRecon_Cleaned.jl`) is publicly available (s.b.). The data for the in-vivo measurement (SPIDI_0007 and SPIDI_0011), run via `RunReconLoop.jl` and `SpiralRecon_SingleIntlv.jl` can be obtained from the authors.

## Getting Started

1. To get started, make sure you have Julia installed. At least v1.6 is preferable (v1.7 is optimal).
2. Clone the GIRFReco.jl project via Github to a local directory
3. Download the data supplement from Zenodo (https://doi.org/10.5281/zenodo.6510020) and extract somwhere. 
   - Note: This might take a few minutes. You can skip ahead and continue with step 5-10 in the meantime.
4. Move the extracted folder (data) into the MRIRecipes.jl directory
5. Clone the `MRIGradients.jl` package by typing `]add https://github.com/BRAIN-TO/MRIGradients.jl.git`
6. Open a Julia REPL in your editor of choice (we recommend VS Code with the Julia extension)
7. In the REPL, type `]` to enter package mode
8. type `activate .` to activate a new Julia environment for the MRIRecipes.jl project
9. type `dev MRIGradients` to tell Julia which MRIGradients to use (this will be fixed upon package registration)
10. type `instantiate` to download and install all of the necessary packages.
11. Proceed to run the demos found in the /recon/ directory.

## Examples
    
1.  **Phantom data** (*presented at [ISMRM 2022, p.2435](https://archive.ismrm.org/2022/2435.html)*): Interplay of GIRF.jl and MRIReco.jl: Using the julia recon package to reconstruct a spiral image from a structure phantom with a GIRF-predicted trajectory 
    ```
    recon/SpiralRecon_Cleaned.jl
    ```   
2.  **In-vivo brain data** (spiral diffusion, *Submitted to ISMRM 2023 and the Sedona Workshop on Data Sampling and Reconstruction*): Interplay of GIRF.jl and MRIReco.jl: Using the julia recon package to reconstruct a spiral image in-vivo (brain diffusion image) with a GIRF-predicted trajectory
    - To run an example image reconstruction (one volume, no diffusion weighting): Execute
    ```
    recon/SpiralRecon_SingleIntlv.jl
    ```
    - To reconstruct all images in a scan, run this loop over all averages and diffusion directions:
    ```
    recon/RunReconLoop.jl
    ```
3.  If you want to try out other example datasets, just rename any of the `ReconConfig_*.jl` files (e.g., `ReconConfig_SPIDI_0011.jl`) to `ReconConfig.jl` and repeat step (2).
