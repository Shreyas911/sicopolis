How to create the Doxygen documentation
---------------------------------------

1. Copy the config template:  
   `cp doxygen1.8.1-config-template.txt my_doxygen1.8.1-config.txt`

2. Edit my_doxygen1.8.1-config.txt:  
   Search for "Revision xxxx", and replace xxxx with the current
   revision number  
   (execute `./rev_id.sh` in the runs directory to find out).  
   Search for "/home/username/Documents/sicopolis/src", and replace
   this by the actual path of the src directory. 

3. Execute the command
   `doxygen my_doxygen1.8.1-config.txt`
   (requires Doxygen version 1.8.1 or later).  
   This creates the html documentation in doxygen/html.  
   Main file: index.html.

4. Documentation as PDF (optional):  
   Go to doxygen/latex and execute `make` (requires LaTeX).  
   This creates the file refman.pdf.
