# phase-field-models

This repo is for developing pf models from scratch and improving understanding from perspectives of mathematics, computer science and material science.

Language: C++, Python
C++ does the compuation and python does the visualization

Reference of numerical Model:

- https://www.youtube.com/watch?v=EZy9lJmMBLs&list=PLaikkcmmwwxL4Ujl4wjPoPp_AHjVZgigc
  The course is made by Prof. M. P. Gururajan (https://www.iitb.ac.in/mems/en/prof-m-p-gururajan), who learns PF from Prof. Peter. W. Voorhees (Northeastern Univeristy) https://pages.nist.gov/pfhub/wiki/voorhees-lectures/

- Programming Phase Field Models, by S. Bulent Biner

- OpenPhase Library, by ICAMS https://openphase.rub.de/index.html

The principles of developing reusable, expandable, managable problems are:

- Do not repeat yourself
  when a block of code are used in different places, put it in a sigle file and call it from other files
- Modularization and Object Orientation Programming
  This is for building new models more effectively and make degugging easier
- ...
