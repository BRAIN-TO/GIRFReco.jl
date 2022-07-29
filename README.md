# MRIRecipes.jl

Recipes for image reconstruction using the MRIReco.jl package. 

Includes a working example for spiral reconstruction with GIRF correction and a Cartesian reconstruction example for sensitivity and off-resonance map calculation.

### NEW DOCS:

Updated to be included in MRIRecipes.jl https://github.com/BRAIN-TO/MRIRecipes.jl

## Getting Started!

1. To get started, make sure you have Julia installed. At least v1.6 is preferable (v1.7 is optimal).
2. Clone the MRIRecipes.jl project via Github to a local directory
3. Download the data supplement from Zenodo (https://doi.org/10.5281/zenodo.6510020) and extract somwhere. 
   - *Note: This might take a few minutes. You can skip ahead and continue with step 5-10 in the meantime.*
4. Move the extracted folder (data) into the MRIRecipes.jl directory
5. Clone the `MRIGradients` project via Github to your Julia dev folder (usually `/home/.julia/dev/` on Linux or `C:\Users\<username>\.julia\dev` on Windows). `MRIGradients` can be found at: https://github.com/BRAIN-TO/MRIGradients
* or you can do this in one step by typing `]add https://github.com/BRAIN-TO/MRIGradients.jl.git`
6. Open a Julia REPL in your editor of choice
7. In the REPL, type `]` to enter package mode
8. type `activate .` to activate a new Julia environment for the MRIRecipes.jl project
9. type `dev MRIGradients` to tell Julia which MRIGradients to use (this will be fixed upon package registration)
10. type `instantiate` to download and install all of the necessary packages.
11. Proceed to run the demos found in the /recon/ directory.

### OLD DOCS

2. Open the GIRFReco folder in the text editor of your choice (we use Atom with the Julia extension *Juno*) after configuring the text editor to use Julia. 

3. Once your editor is prepared, and the GIRFReco folder is opened, open the package cmd line by typing 
   ```
   ]
   ```
   - Your prompt should now say `pkg>` instead of `julia>`
   - All the following commands are within the `pkg>` prompt.
   - If you have to get back to the julia prompt later, press `CTRL+C`
4. Within the package cmd line, activate the project environment:
    ```
    activate .
    ```
    - This activates the Julia environment for the GIRFReco code. 
    - After typing this, wait for a while until the command line finishes setting up the Julia environment. 
5. If this is your first timing activating this environment, you have to instantiate it after activating it via
    ```
    instantiate
    ```
    - Some packages might fail to install, but they are usually not needed by our code.
6. Add MRIReco.jl to your environment's package list via
    ```
    add MRIReco
    ```

## Examples

- To run the following examples, open the listed files in your editor. 
- Make sure your current folder is `GIRFReco` and its environment is activated.
- For Atom/Juno:
     - right-click on the folder, then `Juno -> Work in Folder` and `Juno -> Activate Environment in Folder` (or activate via the `pkg>` prompt as above)
     - Juno shortcut for running an open file: `CTRL+SHIFT+ENTER`
     - Juno shortcut for running a block (separated by a line startin with `## ` in the file): `ALT+ENTER`
    
1.  GIRF.jl in action: Predicting an actual gradient waveform from a nominal one using the GIRF
    ```
    girf/GIRFDemo.jl
    ```
2.  Interplay of GIRF.jl and MRIReco.jl: Using the julia recon package to reconstruct a spiral image with a GIRF-predicted trajectory
    ```
    recon/SpiralRecon.jl
    ```
    
