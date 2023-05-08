# SheffieldAutismBiomarkers

Welcome to the repository of code assocaited with the Sheffield Autism Biomarkers project. 

In the near future the preprint link will be found here: 


And in the somewhat later but hopefully not that much later future the peer reviewed and published link will be found here: 

**CODE:**
This repository contains scripts and functions that fall into several different categories. A detailed explanation of how these different code categories interact with each other can be found in the supplement.pdf file contained in the repository. Here, we provide brief descriptions of the different categories that the code contained in this repository belongs to:  

Core computation functions (layer 4 in the supplement): 
These functions perform direct extraction of key output variables from the EEG data itself. They are: chanEntropyVals.m, chanFuzEnt.m, chanPAC.m, chanPower.m, chanSlopeAlpha.m, getISPC.m, convertCoordinates.m, removeNoiseChansVolt.m, and fileSplitter.m. These files are likely to be the easiest to directly apply to others' data, and therefore they are the most useful out of the box. Each function has comments indicating its necessary input format and what type of outputs it can provide. In some cases, these functions depend on other subfunctions contained in the repository (e.g. chanEntropyVals depends on msentropy.m)

Pipeline function (layer 3 in supplement): 
A pipeline function, as we use the term, is a function that receives input in the form of file paths and metadata directing it to the data that are to be analyzed. It contains calls to the core computation functions. There is one example pipeline function contained in the repository: setReadInAggregate.m. This function will likely not be useful for others to apply to their own data, but it will serve as a useful example for how to structure an analysis pipeline and it demonstrates how data were analyzed for this project. 

Pipeline wrapper script (layer 2 in the supplement): 
A pipeline wrapper, as we use the term, is a script that organizes the metadata for a large-scale analysis and then calls the pipeline function to perform the analysis. There is one example pipeline wrapper script in the repository: getBioConsortDat.m. Note that this wrapper is specific to the Autism Biomarkers Consortium for Clinical Trials dataset. Similar pipeline wrappers were written for other datasets. This script is highly specific to the dataset it was designed to analyze and serves as an example more so than as useful code. 

Bash job initiation script (layer 1 in the supplement): 
The bash job initiation script is a script that asks a high power computing system for resources and submits jobs to be processed in parallel by starting Matlab and calling on the pipeline wrapper script. There is one example bash job initiation script in the repository: exampleBashScript.sh. As with the pipeline wrapper script, this is intended as an example rather than as useful code. 

Auditing code: 
There are two functions which are written to audit progress on a channel by channel phase amplitude coupling analysis: PACaudit.m and audit_singleChanAll.m. 

Summary variable extraction:
The scripts stitchFilesTogether.m and extractingFinalVariables.m were used to combine single channel output variables into participant level output variables and extract final summary variables from participant level outputs, respectively. These are highly specific to the analysis performed for this project and will be most useful as examples rather than immediately useful code. 

Inferential statistics: 
The R script effectSizeCalculations.R was used to perform all inferential statistics and generate many of the plots included in the final paper. It is highly specific to the particular set of analyses performed here. 

**DATA:**
In addition to code, this repository contains the extracted summary variables that were analyzed in this project. Specifically, there are three key .csv files: 
autismBiomarkersAllData2.csv this file contains all of the extracted summary variables with one participant per row of the .csv file. 
testSet.csv and trainSet.csv are the randomly split halves of the full dataset. The code used to perform the randomization can be found in the effectSizeCalculations.R script in a long comment. The randomization was performed only once. 

finally, the file standardEEGlocs.csv contains the standard 32 channel montage that all data were projected into for all of these analyses used throughout this project.
