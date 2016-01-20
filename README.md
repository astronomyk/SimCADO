###Hello world ###

# The README file
###Trying out Git

* Create a file in a folder
* Move to that folder

Tell git about your git account
Please tell me who you are.

<<<<<<< HEAD
```
>> git config --global user.email "you@example.com"
>> git config --global user.name "Your Name"
```
=======
    >> git config --global user.email "you@example.com"
    >> git config --global user.name "Your Name"
>>>>>>> d5630ffe4a9d0f595b050551a2ebbea3068ff960


In the folder

<<<<<<< HEAD
```
>> git init
>> git status
```

Let git know that this file should be tracked


`>> git add readme.md`

Now commit the file. You have to add a comment to the commit

`>> git commit`

To re-upload

```
>> git add readme.md
>> git commit
```
or if I want to recommit all altered files - watch out, dangerous

`>> git commit -a`
=======
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
>>>>>>> d5630ffe4a9d0f595b050551a2ebbea3068ff960


###Branching

List all branches

<<<<<<< HEAD
`>> git branch`

Add a branch

`>> git branch test_branch`
=======
    >> git branch

Add a branch

    >> git branch test_branch
>>>>>>> d5630ffe4a9d0f595b050551a2ebbea3068ff960

Switch to a branch

Now any changes will be registered in test_branch

<<<<<<< HEAD
`>> git checkout test_branch`

Merging two branches

```
>> git checkout master
>> git merge test_branch
```
=======
    >> git checkout test_branch

Merging two branches

    >> git checkout master
    >> git merge test_branch
>>>>>>> d5630ffe4a9d0f595b050551a2ebbea3068ff960


###Playing with Github

To pull down a repository from the web, use:

<<<<<<< HEAD
`>> git clone <url from git>`

=======
    >> git clone <url from git>
>>>>>>> d5630ffe4a9d0f595b050551a2ebbea3068ff960

Github-pages can host a series of webpages - accessable through settings

Add the local repository to a github repository
<<<<<<< HEAD

`>> git remote add origin https://github.com/gastronomyk/SimCADO.git`

Pushing this upto Github

`>> git push -u origin master`

Look in to this a bit more
=======

    >> git remote add origin https://github.com/gastronomyk/SimCADO.git

Pushing this upto Github
>>>>>>> d5630ffe4a9d0f595b050551a2ebbea3068ff960

    >> git push -u origin master
