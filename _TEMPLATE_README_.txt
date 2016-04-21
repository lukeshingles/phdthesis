_README_

This is an RSAA Ph.D. thesis template (though also suitable for an Honours or Masters thesis)   that I used for my 2011 thesis. It is based on the excellent work of Josh Rich, to whom I am indebted. I have added some extensions to Josh's work to make my life easier while writing the thesis, and to add value to the PDF copy that everyone will read (while the bound copies gather dust in the library). Some of these new features are:

* the ability to link objects to SIMBAD by clicking on them (useful when you can't remember that paper that cites your favourite object)

	\object{SIMBAD NAME} --or-- \object{SIMBAD NAME}{Name in the text}

	the latter is useful when the object name is complex, e.g. greek letters

* new bibliography style file, astroads.bst, which places a link to ADS after each paper in the bibliography, and truncates author lists after 8 authors so the SDSS papers don't take up half your bibliography :)

* Tag your data and services to the Virtual Observatory to ensure the provenance and repeatability of your science. For example, show exactly which Hipparcos cone search server and catalogue you used by specifying the IVOA identifier of the service:

	\ivoa{ivo://astronet.ru/cas/hipparcos}

	the resulting link is resolved by the US National Virtual Observatory Directory

As well as these features, the template allows you navigate the PDF (references, citations, contents page) using the LaTeX href package. This is very useful for editing! The default colours are somewhat garish, so please customise to your heart's content inside thesis.tex.

The shell script 'topdf.sh' does the job of compiling pdftex and bibtex for you, so no more manually running latex, bibtex, latex, latex ad infinitum. It also cleans up after itself and leaves all the important files in the 'output' directory.

-------------------

_README_		This file
abstract.tex		Abstract
acknowledgments.tex	Acknowledgments
anu-logo-colour.eps	ANU Logo
anu-logo-colour.pdf	ANU Logo
apj.bst			References in the style of ApJ
astroads.bst		Same, but with links to ADS after each reference
ctable.sty		Good looking tables
dedication.tex		Thesis dedication
disclaimer.tex		Legal disclaimer (required)
macros.tex		Useful symbols, shortcuts etc
references.bib		Example BibTeX database
signature.png		Your signature for the digital copies (optional)
thesis.pdf		Finished product (automatically copied from ./output)
thesis.tex		Main structure, template settings
topdf.sh		Script for compiling pdf from source files (via pdftex)

Each chapter and the introduction, conclusion, appendices etc, have been placed into 
a subdirectory to tidy things up. The 'output' directory is populated by topdf.sh each
time it is run. It contains the various output files produced by pdftex (aux, bbl etc),
including the final PDF file (thesis.pdf), which is copied to the top-level directory.

Any questions, please email me <simon.murphy.nz@gmail.com>

Happy writing!

SJM 21/08/2012


