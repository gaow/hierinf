## Test environments
* local OS X install, R 3.5.0
* Fedora 26, R 3.5.0
* win-builder (devel and release)

## R CMD check results
There were no ERRORs, WARNINGs or NOTEs. 

## R CMD BiocCheck results
There were no ERRORs or WARNINGs. 

There were 6 NOTEs:

* "Consider adding these automatically suggested biocViews: ATACSeq" 

  This biocViews category does not suit our package. 

* "Avoid '<<-' if possible (found in 1 files)"   

  We need the operator <<- to store error messages during parallel computing. 
  See file tryCatch-W-E.R 
  
* "Consider adding a NEWS file, so your package news will be included in 
  Bioconductor release announcements."
  
  See file NEWS.md

* "Consider shorter lines; 188 lines (4%) are > 80 characters long."   
  
  We have lots of error and warning messages which are written on one long line
  since we do not want any line breaks in the error and warning messages. 

* "Consider multiples of 4 spaces for line indents, 1138 lines(26%) are not."   

  We use the default line indents of RStudio. 

* "Cannot determine whether maintainer is subscribed to the bioc-devel mailing 
   list\n(requires admin credentials).  Subscribe
   here:\nhttps://stat.ethz.ch/mailman/listinfo/bioc-devel"
   
  Yes, the maintainer is subscribed to the bioc-devel mailing list. :-)

