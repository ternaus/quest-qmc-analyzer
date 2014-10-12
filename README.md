# README #

This README would normally document whatever steps are necessary to get your application up and running.

### What is this repository for? ###

* This repository contains python program that is used to create plots from the data that [QUEST](quest-qmc.googlecode.com)  package generates.
* At a moment this package is in it's alpha stage. If there will be some interest in it, I would clean the code and release stable version.
### How do I get set up? ###

* Download the latest version of the code to you computer:

```
#!bash

hg clone https://vladimir_iglovikov@bitbucket.org/vladimir_iglovikov/quest-qmc-analyzer
```

* Generate the data using [QUEST](quest-qmc.googlecode.com) package. For now constraint is that name of the output file should start with the lattice name, followed by an underscore, so for the square lattice possible output file names that this analyzing package will be able to work with are: square_123123.out square_asfd123asdsdf.out
* Put this output files into the folder that is named exactly as your lattice. So for the square lattice output files should be found in the "square" folder.
* change variable "folder_with_different_models" in the setting.py file. For example, I keep my outputfiles as:
square lattice => /home/vladimir/work/results/square, honeycomb => /home/vladimir/work/results/honeycomb, thus in the file settings.py I have:

```
#!python

folder_with_different_models = 
```


* How to run tests
* Deployment instructions

### Contribution guidelines ###

* Writing tests
* Code review
* Other guidelines

### Who do I talk to? ###

* Repo owner or admin
* Other community or team contact