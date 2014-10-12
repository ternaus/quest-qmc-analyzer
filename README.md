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

* Generate the data using [QUEST](quest-qmc.googlecode.com) package.
* Put this output files into the folder that is named exactly as your lattice. So for the square lattice output files should be found in the "square" folder.
* change variable "folder_with_different_models" in the setting.py file. For example, I keep my outputfiles as:
square lattice => /home/vladimir/work/results/square, honeycomb => /home/vladimir/work/results/honeycomb, thus in the file settings.py I have:

```
#!python

folder_with_different_models = '/home/vladimir/work/results'
```


* Start the analyzing script:
syntax will be clear from examples:


Examples:

python analyzer.py -m square -u 4 -t 1 -beta 4 -x_variable mu -y_variable rho


This will give you a plot of rho vs mu for square lattice for U =4, t = 1, beta = 4


----


python analyzer.py -m square -u 4 -t 1 -beta 4 -x_variable rho -y_variable m2 -to_screen True


This will give you magnetization squared vs rho for square lattice for U =4, t = 1, beta = 4 and will print array with x_variable, array with y_variable, and array with errorbars to screen

----
python analyzer.py -m square -mu 0 -t 1 -u 4 -x_variable beta -y_variable energy -vline 12

This will plot total energy vs beta for the square lattice for mu = 0, t =1, u = 4 and will plot vertical line at beta = 12
