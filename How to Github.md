References used: 
1. https://docs.github.com/en/get-started/quickstart/hello-world


Github is where we can change and edit the same code files individually while having all of our files tracked and saved online. The purpose is to organize our code and to be able to easily manage the code files without running into problems in the future. Github is a place that stores and organizes our project onto the cloud, so anyone with internet can access it any time. This also means that any changes we want to make requires internet.

To actually use Github as a team, we need to download the Github app onto the desktop or get a Git software.

## What is 'Git' and 'Github'?
Git is the software that tracks the history of changes that people and the team makes on projects. It is different from Github. Github is the internet storage place for our project and provides tools to help teams to collaborate better like command line features, issues (threaded discussions), pull requests, code review, or add-in apps. Github uses Git to track our changes, so that is why we need to install a separate Git software to use Github. In the references used lists Git command lines to use in applications like GitHub Desktop app to get started.

A Github workflow is listed below on the basic workflow when using Github. It starts by creating a repository, and then cloning the repository to your computer, creating a new branch, committing changes, pulling requests, peer reviewing, and then editting the final document.

The first steps are the create and clone a repository.

The actual editting and changing the code part is creating a branch, commit changes, pull requests, and committing to the main document.

## Repository
References used: 
1. https://docs.github.com/en/get-started/using-git/about-git
2. https://docs.github.com/en/get-started/quickstart/github-flow

The starting point of Github is the repository. The repository is the location of all our of project files. It also shows all of the changes and updates made to each file. Right now- 'thermal-matlab' is our repository name and where all the documentation and files are in. Repositories can contain folders and files, images, videos, spreadsheets, and data sets -- anything your project needs. The respository is already made for this project but for reference in the future, below is how you make a repository.

Every repository has a README.md file. The README.md file is to outlines the information of our project and code files. 

To create a repository:
1. In the upper-right corner of any page, use the  drop-down menu, and select New repository
2. In the Repository name box, enter a name.
3. In the Description box, write a short description.
4. Select Add a README file.
5. Click Create repository.

## Cloning the Repository
References used: 
1. https://swcarpentry.github.io/git-novice/08-collab/index.html
2. https://docs.github.com/en/repositories/creating-and-managing-repositories/cloning-a-repository

This means to copy the repository onto your computer locally so you can make changes, edits, add and remove files, and push larger commits. When you clone the repository, you download the full copy of the repository at that point in time. If other people made and committed changes, you can add those changes to your copy of the repository. 

Changes to the text files like this one and the README doesn't require cloning, and can be editted here. But with files that use a different software, in our case, MATLAB, we may want to consider making a copy of this repository onto your computer. This requires to download an extra app onto your computer. 

This is where 'Git' itself comes in and requires each person to have their own verison of Git downloaded. Github reccomends using their Github Desktop App, but you can download a Git software like Gitbash instead. I'm going to only go over how you can clone your repository on Github Desktop App, but there's other references if you want to use something else.

## Github Desktop app, Configurations and How to Clone the Repository
References Used: https://docs.github.com/en/desktop/installing-and-configuring-github-desktop
Install Github Desktop App with this guide here, just follow how to install it and adding in your account: https://docs.github.com/en/desktop/installing-and-configuring-github-desktop/installing-and-authenticating-to-github-desktop/authenticating-to-github

With Github Desktop App installed, this is how you clone the repository:
1. Go to Files -> Clone Repository
2. Edit to where you want it stored on your computer
3. Save, should be a file on your computer :) 
4. If you don't see the repository, make sure you're using the correct account that the repository is linked to.

### A few features on the Github App: (all the weird terminology is explain below)
#### "Pull Origin"
Its the rightmost tab on the app's interface. As commits or changes are done on the branch, you can keep your local copy of the project in sync my pulling from the remote repository. As you sync, the branch that you're in updates. If you want to change branches, click the drop down menu under 'current branch' and change there.


## Branches 
The main purpose is to avoid saving different verisons of the same file and branches accomplish similar goals in GitHub repositories. Branches also enables your collaborators to see your ongoing work too.
Something like this can be easily avoided by using branches:

story.txt

story-joe-edit.txt

story-joe-edit-reviewed.txt

Github uses branches for keeping bug fixes and feature work separate from our main (production) branch. When a change is ready, they merge their branch into the 'main' branch. When you make a new branch, you can do whatever you want with the repository (add and edit files, rename a file, move files, delete anything). In a new branch will not change anything in the main branch.
Everyone starts off working in the main branch which is where all edits are officially published.

When you create a branch off the main branch, you're making a copy, or snapshot, of main as it was at that point in time. If someone else made changes to the main branch while you were working on your branch, you could pull in those updates. Everytime you want to change the code file or work on it, create a branch. If everyone wants to work on the same file and makes their own branches, it won't affect the file at all and as many people can work on it at once. Once the workflow is complete, the branch can be deleted.

To create a branch:
1. Click the code tab in repository
2. Click the drop down menu at the top of file that says main
3. Type in a branch name into text box
4. Click create branch

## Committing and Pushing Changes
This is where all the edits happen. 

In a new branch (where it is just you making edits thats not affecting what is already published or being editted at the same time) you can make changes to the files in your repository. Every time you edit and want to make a final change to the file, Github is set up that after you finish editting, you can write a commit message, describing everything you changed and why. They help capture the history of the what changes and who changed what and why it was changed.

Commit means to save your edits into your local repository. Ideally, each commit contains an isolated, complete change like 'fix typo' or 'increase rate limit'. This makes it easy to revert your changes if you decide to take a different approach. For example, if you want to rename a variable and add some tests, put the variable rename in one commit and the tests in another commit. Later, if you want to keep the tests but revert the variable rename, you can revert the specific commit that contained the variable rename. If you put the variable rename and tests in the same commit or spread the variable rename across multiple commits, you would spend more effort reverting your changes.

Push means to transfer your commits to the remote location of reporsitory. This means that you can access your work from any device. It also means that your collaborators can see your work, answer questions, and make suggestions or contributions. 

Continue to commit and push till ready to peer review.

To commit a change in the GitHub Webiste:
1. Open the file and click edit (pencil icon)
2. Edit file
3. At the bottom of page, there is a commit changes box to add comments
4. Click commit changes

To commit and push a change on your computer (using Desktop App):
1. From your cloned repository, open the file
2. Edit File normally
3. Save and go back to desktop app and add comments and click 'commit to 'branch name''
4. Click push changes too

## Pull Request
Referenced used: 
1. https://docs.github.com/en/pull-requests/collaborating-with-pull-requests/proposing-changes-to-your-work-with-pull-requests/creating-a-pull-request
2. https://docs.github.com/en/get-started/writing-on-github/getting-started-with-writing-and-formatting-on-github/basic-writing-and-formatting-syntax
3. https://docs.github.com/en/pull-requests/collaborating-with-pull-requests/collaborating-on-repositories-with-code-quality-features/about-status-checks

When we are done comitting, its time to peer review. Creating a pull request just means to ask collaborators for feedback on your changes. Pull request reviews are considered the most important part of the process and can be so valuable that some repositories require an approving review before pull requests can be merged. Pull Request can be done on the website or the Github app (not applicable to our project).

Steps to make a pull request:
On GitHub.com, navigate to the main page of the repository.

1. In the "Branch" menu, choose the branch that contains your commits.
2. Above the list of files, click  Pull request.
3. Use the base branch dropdown menu to select the branch you'd like to merge your changes into, then use the compare branch drop-down menu to choose the topic branch you made your changes in.
4. Type a title and description for your pull request.
5. To create a pull request that is ready for review, click Create Pull Request. To create a draft pull request, use the drop-down and select Create Draft Pull Request, then click Draft Pull Request. For more information about draft pull requests, see "About pull requests."

On the Desktop App:
1. Click Create Pull Request. GitHub Desktop will open your default browser to take you to GitHub.
2. 


When you create a pull request, include a summary of the changes and what problem they solve. You can include images, links, and tables to help convey this information. If your pull request addresses an issue, link the issue so that issue stakeholders are aware of the pull request and vice versa. Also in your pull request, you can add comments to specific lines of the code and @ people who you want to see this review. 

See Reference 2: Basic syntax of pull request for your own knowledge, but not used for this project
See Reference 3: you can configure your repository to have a check list on pull requests, you will see any checks that failed on your pull request. This helps you catch errors before merging your branch.

## Finishing up a round of editting
After creating a new branch, committing/pushing changes, and finishiing/merging a pull request, then the cycle of editting a file is completed and here are the final steps listed.
1. delete your branch after the pull request is completely finish. This prevents anyone else messing with it by accident when everything is already done. The history from that branch is already recorded and track so no information is lost.

Once this branch is deleted, time to complete a whole new cycle of edits.






Reference used to edit .md files: https://docs.github.com/en/get-started/writing-on-github/getting-started-with-writing-and-formatting-on-github/basic-writing-and-formatting-syntax
