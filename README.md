This is a collection of tools meant to be used for analysis of disperse multiphase results from Bluebottle simulations.
All of the useful ones are written in c and Python, everything else is deprecated.
Feel free to use and modify these under the terms of the Apache License included in this folder.
The most useful way to use the tools is to probably create a fork of this repository and work on it there.
If there are any bugs or additions you feel would be useful to any users, please create a pull request.

Note that all of the tools have been developed for use with triply-periodic sedimentation (gravity is positive in -z) -- this means that I don't have explicit checks anywhere for boundary conditions!!
Different simulations will probably require some modifications for boundary conditions, if you want to make the simulations general enough for this please feel free and then submit a pull request!

I apologize in advance for some of the quirks -- none of this was developed with a broad "plan" in mind, nor was it developed to be production code (as opposed to the main bluebottle repository).
As such, there are many areas that require file structures to be a certain way, files to be in a certain place, etc.
I've tried to write down as many of them as I remembered, but I'm sure I missed several.

On a similar vein, this repository is CONSTANTLY changing (including file structure updates, I seem to love doing those every month or so...), so check back often for bug fixes, extra features, free ice cream, and other goodies.

# Main Directory
There are several top-level directories. 
I use one of them.
* bbtools
  * c-tools
    * This is where every useful tools resides and all of the current development is taking place
  * deprecated
    * Stuff that I've moved around but not wanted to throw out yet that I either don't use or may come back to
  * matlab
    * OLD matlab tools. I think every one of these has been replaced by a c tool. MATLAB is slow. I don't use MATLAB anymore.
  * presentation
    * ignore this
  * sppressions
    * suppresion file for using with valgrind when checking out cuda stuff

# File Structure
The tools are meant to be used on a standard Bluebottle simulation directory with this structure:
* Top-Level-Simulation-Directory
  * input
    * flow.config
    * part.config
    * record.config
  * output
    * flow-%lf.cgns
    * part-%lf.cgns
    * grid.cgns
  * record
    * Files are not used in analysis

In addition, the following analysis directory has been added to hold the results of the analyses:
* Top-Level-Structure-Directory
  * analysis
    * fourier-reconstruction
    * part-pair-distribution
    * phase-averaged-velocity
    * tetrads

I created this structure because it was convenient at the time, so all of the analysis relies on the structure to run correctly.
If you'd like to take it upon yourself to make it more general, that'd be great -- but submit a pull request!
This structure also mimics (I think -- there may be one or two differences that I've forgotten to change) the structure that this repository has taken.

# Tools
The following tools are (currently, as of 08/03/2016) included:
* fourier-reconstruction
  * 1-dim-flow
    * Reconstructs a fluid field averaged out over the cross plane (x,y) directions. Uses part-%lf.cgns files (essentially reconstructs using delta functions where particles are)
  * 1-dim-part
    * Reconstructs a particle psuedo-field averaged out over the cross (x,y) directions. Uses part-%lf.cgns files (essentially reconstructs using delta functions where particles are)
  * 3-dim-part
    * Reconstructs a particle field in 3-dim using the phase array from flow-%lf.cgns. It would be possible to use this to reconstruct in 1-dimension as well. (essentially reconstructs using Heaviside functions where the particles are)
* part-pair-distribution
  * 3d-hist-bins
    * Conditional pair distribution function based on volume fraction (bins)
  * part-pair 
    * Several different methods, including binning, for calculating the pair distribution function. Right now only binning seems to converge nicely, unfortunately.
* phase-averaged-velocity
  * flow
    * Phase-averaged flow velocities over time
  * part
    * Part-averaged flow velocities over time
* tetrads

There are several other subtools present in some of the directories.    

# Tool Directory
Each tool directory is based around the following structure:
* Top-Level-Tool-Directory
  * Makefile
    * self-explanatory
  * src
    * all source code.
  * sim
    * sim binary and config file
  * plt
    * python plotting tools

Running a simulation is usually done from the sim directory using "./binary /path/to/the/top/level/simulation/directory".
It is necessary that within the simulation directory, in the correct folder, a config file exists. 
The correct folder is the one given by, for example: "..../analysis/fourier-reconstruction/3-dim-part/config", although the name
of the config file there is not correct.
All config files are initially located in the tool sim directory, simply copy them where you want them.
If not, the program will fail.

I'm sure there are several other quirky things going on, but there's a fair amount of error checking, so if files aren't found the program should fail.
