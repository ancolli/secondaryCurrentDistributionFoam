# v1.0.0
# secondaryCurrentDistributionFoam
It is described how to simulate secondary current distribution in electrochemical reactors with a multiregion approach (solve solid and fluid phase) whit the help of OpenFOAM, the developed solver multiregionPotentialFoam and the BC regionCoupleSolidFluid. It shows how to pre-process, run and post-process a basic case involving in a 2D domain. 
The proposed strategy allows to have as an input i-th electrochemical reactions per electrode, including bipolar electrodes in which at least two reactions take place (reduction and oxidation).

# Disclaimer
This offering is not approved or endorsed by OpeFOAM Foundation, producer and distributor of the OpenFOAM software via www.openfoam.org.

# Usage
In applications (A) you will find the scripts to compile the solver, BC and post-processing utility in order to solve secondary current distribution.
In tutorial (B) you will find an example of a 2D cell composed of 1 cathode, 1 anode and 1 bipolar electrode. 
In more information (C) you will find a brief description of the proposed tool

# #  A) Applications
1)  	
A- Paste applications/utilities/Solvers/multiRegionPotentialFoam inside OpenFOAM user directory (Applications/Utilities/Solvers).  
B- Open a terminal inside multiRegionPotentialFoam.  
C- Run wmake.  
2)	
A- Paste applications/utilities/BC/regionCoupleSolidFluid inside OpenFOAM user directory (Applications/Utilities/BC).  
B- Open a terminal inside regionCoupleSolidFluid.  
C- Run wmake.   
3)	
A- Paste applications/utilities/postProcessing/potentialWallFlux inside OpenFOAM user directory (Applications/Utilities/postProcessing).
B- Open a terminal inside potentialWallFlux.
C- Run wmake.   


# #  B) Tutorial
1- Paste tutorial inside OpenFOAM user directory (Run/Tutorials).
2- Enter to tutorial and open a Terminal.
3- Modify properties inside system/anode-cathode-bipolar-electrolyte/changeDictionaryDict.
4- Run ./Allrun.


# # C) More information
For a multiregion partitioned solver the working steps are as follows:
1. Define multiple meshes, one for each "region".
2. Create field variables on each mesh.
3. Solve separate governing equations on each mesh.
4. Multiregion coupling at the boundary interface between regions.
5. Subiterate until fully-coupled solution is reached.
Case Setup
In this case we have four different regions Liquid1 (electrolyte), Solid1 (cathode) Solid2 (bipolar electrode) and Solid3 (anode). The case is a simple set up with fixed different potentials for the terminal electrodes (solid parts) (Solid1 and Solid3), which are surrounded by the electrolyte (Fluid1).

Geometry and Mesh
The geometry is defined and then meshed using the OpenFOAM blockMesh tool. After running the blockMesh utility. The regions are created in the domain depending upon their phases and they are createdon the basis of the zones defined.
In the presented case the following regions must be created: electrolyte, cathode, bipolar and anode. In order to create these regions the domain is divided into zones. To do so, a subset of cells within the domain is selected to form a so-called cellSet. In order to define cellSets and cellZones a OpenFOAM command line utility called topoSet is used. This tool requires a dictionary-file as input. Thus, within the case folder a file topoSetDict must exist, which contains the settings used to define the different cellSets/cellZones/regions in the domain.

Splitting the Mesh
After the user has defined all necessary regions by creating zones for them as described in the previous section the mesh of the domain has to be split into several disjoint meshes. Note that the originally created mesh of the full domain will be used within the regions.

Files setup by the user
When starting a new multiregion case the directories and their content must be created manually by the user according to the problem definition.

0 directory: First, manually bring in the necessary field files as usual. For a multiregionPotentialFoam case it is necessary to have files for: fi and k for each phase. 

constant directory: As is any standard OpenFOAM case, the constant folder must contain a standard polymesh directory, generated by the standard blockMeshDict-file, which defines the full domain and its mesh and is located in the system directory.

In contrast to a standard case, the files defining the other properties have to go into the different regional folders, i.e. transportProperties within the folders for fluid regions.
Within the constant-folder it is necessary to produce all the region folders. Additionally within the constant folder it also is necessary to build (or copy from elsewhere) a file called regionProperties. This file assigns the physical phase to each region: Either fluid or solid.

system directory: In the system directory, once again set up the folders for the regions. In each region folder there should be a changeDictionaryDict file, which contains details about the necessary fields in the region like fi, and k.

Get a working controlDict file, for example from the tutorial into this folder. Afterwards get a dummy fvSchemes file. This one is the same as any other fvSchemes, except for the different functions containing no values between the curly brackets. Further on one has to get a fvSolution which only defines the outer correctors into this folder. It is optional to get a decomposeParDict file in case one opts for running parallel computations (not included in the example). For all of the different regional folders: Get a decomposeParDict file and get full fvSchemes and fvSolution files into the folders. For the latter ones, keep in mind that they will be different for the fluids and for the solids.

Files setup by OpenFOAM utilities
Most of the necessary case files and folder for the different regions are created by automated generation using scripts and OpenFOAM utilities. The most important OpenFOAM utilities for a multiregion case are:
splitMesh creates the polyMesh directories and their content within the constant/region_Solidi-Fluid folders; additionally it creates 0/region_Solidi-Fluid/ directories for all regions and copies all the field files existing in the 0 directory into the 0/region_solid-fluid/directories
changeDictionary uses changeDictionaryDict files located in system/region_Solidi-Fluid/folders to create initial, boundary and coupling conditions for all fields existing in 0/region_Solidi-Fluid/ directory for all regions

Post-processing
By the command line >> multiRegionPotentialFoam -postProcess -func surfaces -latestTime it is possible to obtain in the region: electrolyte, and by means of the post processing utility potentialWallFlux (supplied here) the total current and current distribution per electrode.

Running the case
After following each step as defined in the above tutorial. The case can be executed by using the command in the terminal ./Allrun. An Allrun can be explained as a script file which contains all the commands used to execute the case.

