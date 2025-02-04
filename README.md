Welcome to my this repo which will be used for Tutorial purposes!

This repository exists to learn how to clone, add files, check status, create branches, commit changes, and push  those changes.

These files are experimental data from Li et al. 2021 (JFM) titled ``Experimental study of a turbulent boundary layer with a rough-to-smooth change in surface conditions at high Reynolds numbers''. This data is from an experiment that investigated internal boundary layers over smooth-to-rough surfaces, controlling for Reynolds number (Re) and roughness (ks). 

Two sets of datasets are included in the ``ExperimentalData'' directory: ReXXks16, where the sandgrain roughness height of the downstream roughness is held constant in order to investigate the effects of Reynolds number on the downstream adjustment to equilibrium, and ks_Re14ksXX where the Re is held constant and the sandgrain roughness height is  adjusted. 

These directories include both the original data in .txt form, as well as the cleaned-up data to be used in with the script (i.e., Re07ks16_xhatXX, where xhat is the downstream distance from the roughness transition). Additionally, we have files labeled *_BL which have characteristics of the upstream boundary layer.

												     
#					COMMAND LINE TUTORIAL					     #
												     


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
This will add changes to all files in the directory

	git add . 
 
OR this only adds the changes made to this file

	 git add Li_VelScaling.m 

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


#					GitHub Desktop Tutorial					     #


So, you have an aversion to command lines? No matter, you are in luck! GitHub has a desktop app that makes all of thisnonsense typing go away, and replaces it with simple point and clicks!

First, you'll want to download GitHub Desktop:

	https://desktop.github.com/download/

Next, you will sign-in and see a screen like this (yours will be a little different if this is your first sign-in).

![Alt text](./Images/StartingScreen.jpg?raw=true "Start")

Click the dropdown (shown with the red box) to bring up the add -> clone repository option

![Dropdown](./Images/AddRepo.jpg?raw=true "Add")

Next, click the clone repository button and you will see the box below, click url and add the following URL there:

	https://github.com/JCooke188/MyTutorial.git

Make sure you select where you want to create the cloned repository!

![Clone](./Images/CloneRepo.jpg?raw=true "Clone")

Once you've hit 'Clone' you will see something like this:

![Inside](./Images/IntheRepo.jpg?raw=true "Inside")

Find the files in the finder (or if on Windows the File Explorer) and make a change to the Li_VelScaling.m matlab script like:

	% Edited by: (your name) on (today's date)

Once you've done that, save your file and go back to GitHub Desktop.

![Commit](./Images/CommitChange.jpg?raw=true "Commit")

You will see the name of the file we changed on the left side, as well as the reflected changes on the right. In the bottom lefthand side, create a commit title like 

	Added the date

As well as a note 

	Today is Feb 4th.

Then 'Commit to Main'. After you've clicked this, you should see this screen:

![Push](./Images/PushRepo.jpg?raw=true "Push")

Push with the command+P or ctrl+P shortcut, or simply click "Push origin". Congrats, you are now an expert at using GitHub Desktop!

Thank you for checking out this tutorial! Best of luck on your git journey.
