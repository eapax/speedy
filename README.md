# speedy
A branch of SPEEDY ver 41, forked from [Adam Paxton](https://github.com/eapax/speedy).

Please see the original repo and references therein. 

Changes in this branch compared to the original can be see at the [git compare](https://github.com/eapax/speedy/compare/stochastic_rounding...tomkimpson:stochastic_rounding)

There is also a draft PR here https://github.com/eapax/speedy/pull/1/files

note we have filtered to just show changes to source code (`.f90`) - does not show any changes to e.g. namelist files.


Major changes

geop

removed most set precision options - now set once globally, with the exception of the following:



Minor changes

* Remove nc from gitignore so we have copies of the input files in the remote

* Some general deletions in junk etc