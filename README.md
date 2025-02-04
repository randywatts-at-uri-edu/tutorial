Welcome to my this repo which will be used for Tutorial purposes!

This repository exists to learn how to clone, add files, check status, create branches, commit changes, and push  those changes.

These files are experimental data from Li et al. 2021 (JFM) titled ``Experimental study of a turbulent boundary layer with a rough-to-smooth change in surface conditions at high Reynolds numbers''. This data is from an experiment that investigated internal boundary layers over smooth-to-rough surfaces, controlling for Reynolds number (Re) and roughness (ks). 

Two sets of datasets are included in the ``ExperimentalData'' directory: ReXXks16, where the sandgrain roughness height of the downstream roughness is held constant in order to investigate the effects of Reynolds number on the downstream adjustment to equilibrium, and ks_Re14ksXX where the Re is held constant and the sandgrain roughness height is  adjusted. 

These directories include both the original data in .txt form, as well as the cleaned-up data to be used in with the script (i.e., Re07ks16_xhatXX, where xhat is the downstream distance from the roughness transition). Additionally, we have files labeled *_BL which have characteristics of the upstream boundary layer.

######################################################################################################
#												     #
#					COMMAND LINE TUTORIAL					     #
#												     #
######################################################################################################


Clone the repository to your local machine:

	git clone https://github.com/JCooke188/MyTutorial.git

Enter the Repository
	
	cd MyTutorial

Open the README.md (this file!)

	vim README.md

Let's see if we accidentally changed anything

	git status

We should see something like 
	On branch main
	Your branch is up to date with 'origin/main'.
	
	nothing to commit, working tree clean

Now, let's change something. Open the Matlab script -- you can use vim, VSCode, Matlab, etc. since we are opening it to change it. Under '%% Start' let's add:

	% Edited by: (your name) on (today's date)

Once that change has been made and saved on the file, try git status again, we should now see that instead of 'nothing to commit, working tree clean' it says:
	Changes not staged for commit:
	(use "git add <file>,,," to update what will be committed)
	(use "git restore <file>..." to discard changes in working directory)
		modified:    Li_VelScaling.m

Now, we will add our changes to be staged for committing.

	git add . (this will add changes to all files in the directory) OR
	git add Li_VelScaling.m (this only adds the changes made to this file)

If we do git status again, we will see now
	Changes to be committed:
		(use "git restore --staged <file>,,," to unstage)
			modified:    Li_VelScaling.m (in green!)

Let's commit our changed files now:
	git commit -m "Committment isn't that hard with git"

The -m command allows us to add a message to our commit so we know what we committed at the time.

Lastly, we need to push our changes

	git push -u origin main

where we push our branch (main) with -u to add upstream tracking reference to the repository (origin)

 
