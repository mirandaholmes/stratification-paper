# stratification-paper
Code, outputs, data, and figure-making files for "Simulating sticky particles: a Monte Carlo method to sample a stratification"

This folder contains all* the codes required to regenerate the data in the paper, 
as well as the data and program output. Here is a summary of what it contains: 

* Exception: PolymerStrat_s1_labels.txt, PolymerStrat_s2_data.txt, trimer_d2_k1_data2.txt have
not been uploaded yet, because they are too big to upload using the web interface. 


------------------------------------------------ 
Folder CodesOutputs/ 
------------------------------------------------ 

This contains one folder per numerical experiment. Each folder is self-contained, 
and will run on its own (with the necessary libraries.) 
For each experiment, all the parameters are set in file testEXPTNAME.cpp. 

To run a given numerical experiment, 
you must first download the Eigen library, available here: 

http://eigen.tuxfamily.org

Then you must modify the makefile so that the directory given in OPTFLAGS points to 
the location of your Eigen library. After this, each numerical experiment can be run 
by changing to its directory and typing “make”. 

Some experiments have set seeds, in which case you will reproduce the data in the 
paper (assuming the same kind of implementation of the random number generator 
on your system.) Others have the seed set to 0, which chooses a seed based on 
the clock time. 

For the trimer experiment, it was run using different values of kappa. Only the script 
for kappa=1 is included, however you can easily run for other values of kappa by 
opening "testTrimer.cpp" and changing the line "double kappa = 1;" to whatever value 
you desire. 

Also note that for the trimer experiment, we have also included two versions: one which 
lists the set of manifolds in advance (SYSTEM=Trimer in the makefile), and one which 
doesn't list them in advance, but rather calculates them on the fly (SYSTEM=Trimer0 
in the makefile.) We ran the first version to generate data for the paper, but both generate equivalent data and so they are both included as an example. 

All output from each experiment for the data referenced in the paper is in files 
beginning with “output_EXPT.txt”, where EXPT references the name of the experiment.
This output lists the parameters used, as well as the value of the random seed, 
rejection statistics broken up by type of rejection, and the time it took the program to run. 



------------------------------------------------ 
Data/ 
------------------------------------------------ 

This folder contains the data used to make the figures and statistics reported 
in the paper. It was generated using the codes given in CodesOutputs/. 

For the trimer experiment, file "trimer_d2_k1_data2.txt" was run using 10x more points
than the script in CodesOutput/trimer uses, 10^7 points instead of 10^6. 



------------------------------------------------ 
Present folder
------------------------------------------------ 

Contains the Matlab scripts used to make certain figures, as well as the saved 
figures (in epsc format.) 


------------------------------------------------ 
List of Experiments
------------------------------------------------ 

b1-b5: Brownian dynamics experiments, used to compare with data from Example 3. 
Used in makefigs_polymer.m. 


BadMoves: implements the "poor" choices of proposals described in Section 3.3.1. 
No figure; rejections statistics referenced in text. 


e1-e3: ellipsoid sampler, Example 5. 
e1 is a 3d ellipse, no longer referenced in text. 
e2 is a 10-dimensional ellipse, stratification is constructed by slicing with planes. 
e3 is a 10d ellipse, stratification is surface = interior. 
Data analyzed in hists_ellipsoid.m. 


parline: Example 1. Plotted in makefigs_parline.m 


s1-s3: Example 3. 
s1 uses k=2, for extrapolating data. 
s2 uses k=2.885, for comparing efficiency to Brownian dynamics simulations. 
Data from both is plotted in makefigs_polymer.m, and also makefigs_kAB.m (s1, 
to address self-assembly question.) 
s3 has small number of samples; used to make images of polymer, 
in makefigs_polymerimage.m. 


trimer: Example 2. Called with several values of kappa, and # of points, 
as described in text. Figures made in makefigs_trimer.m. 


w1-w22: Example 4, a flexible polymer. Experiments with various parameters. 
Figures made in makefigs_polymerwall.m. 


wb1-wb7: Example 4, semiflexible polymer (has bending stiffness). 
Figures made in makefigs_polymerwall.m, makefigs_polymerABABimage.m. 


wXplot: same as above two items, but run with small numbers of points simply 
to make plots. Plots make in makefigs_polywall_image.m. 


------------------------------------------------ 
------------------------------------------------ 

