# The README file
Trying out Git
*Create a file in a folder

Tell git about your git account
*** Please tell me who you are.

    >> git config --global user.email "you@example.com"
    >> git config --global user.name "Your Name"


In the folder

    >> git init
    >> git status

Let git know that this file should be tracked

    >> git add readme.md

Now commit the file
You have to add a comment to the commit

    >> git commit

To re-upload

    >> git add readme.md
    >> git commit
or if I want to recommit all altered files - watch out, dangerous

    >> git commit -a


Branching
List all branches

    >> git branch

Add a branch

    >> git branch test_branch

Switch to a branch
Now any changes will be registered in test_branch

    >> git checkout test_branch

Merging two branches

    >> git checkout master
    >> git merge test_branch



Playing with Github
To pull down a repository from the web, use:

    >> git clone <url from git>

Github-pages can host a series of webpages - accessable through settings

Add the local repository to a github repository

    >> git remote add origin https://github.com/gastronomyk/SimCADO.git

Pushing this upto Github

    >> git push -u origin master
