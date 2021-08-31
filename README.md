# GIRFReco

Minimum working example for spiral reconstruction with GIRF correction

## Getting Started!

To get started, make sure you have Julia installed. At least v1.5 is preferable (v1.6 is optimal).

Then, open the GIRFReco folder in the text editor of your choice (we use Atom) after configuring the text editor to use Julia. 

Once your editor is prepared, and the GIRFReco folder is opened, open the package cmd line (by typing "]") and type the following:

"activate ."

This activates the Julia environment for the GIRFReco code. 

After typing this, wait for a while until the command line finishes setting up the Julia environment. If there is an error, you might have to instantiate the environment (follow the prompts in the command line)

After instantiation (successful or not), you have to add MRIReco.jl to the package list

## Examples
1.  GIRF.jl in action: Predicting an actual gradient waveform from a nominal one using the GIRF
    ```
    girf/GIRFDemo.jl
    ```
2.  Interplay of GIRF.jl and MRIReco.jl: Using the julia recon package to reconstruct a spiral image with a GIRF-predicted trajectory
    ```
    recon/julia_recon_spiral_w_girf.jl
    ```

