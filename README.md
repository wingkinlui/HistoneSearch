# HistoneSearch
A python package (written in Python 3.8.5) to identify and quantify histone proteoforms, including co-isolated isobaric proteoforms, from top-down MS/MS spectra.
Require numpy, pandas, pickle, pulp, os, scipy and xml

The HistoneSearch tool contains four modules: The major module “HistoneSearch.py”, the spectra file processing module “Deconvolution_File_Processor.py”, the search module “Search_Function.py” and the quantification module “FIRR_Quantification.py”.
The major module imported the deconvoluted spectra file (from either TopFD (_ms2.msalign) or from cRAWler (.puf)), the candidate sequence file (.fasta) and the three other modules. 

The deconvoluted spectra would be loaded to the processing module. All the necessary spectral information would be extracted and saved to a python dictionary. The spectra dictionary would then be loaded to the search module and PrSMs will be returned to each MS/MS spectrum. The PrSMs would subsequently be loaded to the quantification module. Linear optimization would be performed for each MS/MS spectrum for co-isolated isobaric proteoform quantification. Eventually, a .csv table will be generated to report the identified and quantified PrSMs.
