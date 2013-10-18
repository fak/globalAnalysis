# _Global Analysis of Small Molecule Binding to Related Protein Targets_

_Description: Set of scripts accompanying chapter 4 in my PhD thesis: 'Integration and analysis of protein evolutionary relationships and small molecule bioactivity data'; These scripts slightly modified from a [study](http://www.ploscompbiol.org/article/info%3Adoi%2F10.1371%2Fjournal.pcbi.1002333) published in PLoS CB_

## Project Setup

Scripts are written in python - scripts are used to access the database and process the query results. A second layer is written in R for statistical analysis and data plotting. The analysis layer is called from a knitr document in the wiki section.  Several questions are addressed in the study and each has an associated python script and knitr document for analysis:  

1. interAssay.py
2. paralogs.py
3. orthologs.py

## interAssay.py
This script runs through the routines necessary to assess the potency difference observed when comparing measurements taken against identical targets but reported from different assays. 

## paralogs.py
This script runs through the routines necessary to assess the potency difference observed when comparing measurements for the same compound taken against a human target and one of its paralogs.

## orthologs.py
This script runs through the routines necessary to assess the potency difference observed when comparing measurements for the same compound taken against a huma
n target and its rat ortholog.

An overview of the analyses conducted for this study can be found in the wiki section. All R code used to carry out the analyses is specified in the wiki pages.Input data files can be downloaded [here](http://www.ebi.ac.uk/~fkrueger/repo/globalAnalysis/data/). 

## Get in touch

- corresponding author: Dr John P Overington, jpo@ebi.ac.uk
- code and analysis: Felix A Kruger fkrueger@ebi.ac.uk
- twitter: [@mo_sander](https://twitter.com/mo_sander)

## License

Copyright [2012] [EMBL-EBI]

   Licensed under the Apache License, Version 2.0 (the "License");
   you may not use this file except in compliance with the License.
   You may obtain a copy of the License at

       http://www.apache.org/licenses/LICENSE-2.0

   Unless required by applicable law or agreed to in writing, software
   distributed under the License is distributed on an "AS IS" BASIS,
   WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
   See the License for the specific language governing permissions and
   limitations under the License.
