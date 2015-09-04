# zprime13-treemaker

The treemaker for the 13TeV version of the zprime analysis (both channels). 

Probably will be ported to https://github.com/TC01/Treemaker as plugins
later once Ben gets around to it.

This is intended to run over Janos's B2G ttrees, which are almost-direct
conversions of the new 13 TeV ntuples.

# How to use

The ```run.py``` script runs all the right files and outputs trees.
You can modify that to change what files get ran over.

The ```mp_run.py``` is a version of ```run``` that uses multiprocessing,
as an experiment.

If you don't want to the fancy stuff done in either run script, you
can just copy the ```run()``` function into a script of your own and
hard-code what you want to run over.
