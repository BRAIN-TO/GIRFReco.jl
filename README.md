# GIRFReco

Minimum working example for spiral reconstruction with GIRF correction

## Getting Started!

1. To get started, make sure you have Julia installed. At least v1.5 is preferable (v1.6 is optimal).

2. Open the GIRFReco folder in the text editor of your choice (we use Atom with the Julia extension *Juno*) after configuring the text editor to use Julia. 

3. Once your editor is prepared, and the GIRFReco folder is opened, open the package cmd line by typing 
   ```
   ]
   ```
   - Your prompt should now say `pkg>` instead of `julia>`
   - All the following commands are within the `pkg>` prompt.
   - If you have to get back to the julia prompt later, press `CTRL+C`
5. Within the package cmd line, activate the project environment:
    ```
    activate .
    ```
    - This activates the Julia environment for the GIRFReco code. 
    - After typing this, wait for a while until the command line finishes setting up the Julia environment. 
4. If this is your first timing activating this environment, you have to instantiate it after activating it via
    ```
    instantiate
    ```
    - Some packages might fail to install, but they are usually not needed by our code.
5. Add MRIReco.jl to your environment's package list via
    ```
    add MRIReco
    ```
6. We have to make some small changes to MRIReco to make it compatible with GIRFReco. To make packages editable add them as development packages via
    ```
    dev MRIReco
    ```
    - *Note*: If you have installed MRIReco as a dev package for other projects earlier, it might not have been updated to the latest version. In this case, go to you local package folder (usually `C:\Users\<username>\.julia\dev\MRIReco`) and `git pull` to update to the current master
7. Update your packages and restart Julia. You might have to repeat sthis several times until the update command finishes without errors:
   ```
   update
   ```
   - *Note*: The `dev` package of MRIReco might have a local path from a different computer saved in the `Manifest.toml` file (main folder of `GIRFReco`). If you get a related error message, use a text editor to change it in that file to `C:\Users\<username>\.julia\dev\MRIReco` and redo the update.

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
    recon/julia_recon_spiral_w_girf.jl
    ```

