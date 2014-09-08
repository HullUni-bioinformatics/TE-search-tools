#README

# TE discovery in a genome assembly

## Overview

These are a set of wrappers compossing a pipline which finds transposable elements in a genome assembly. The pipline includes the following steps:  

1. **MAKE A DENOVO LIB WITH REPEAT-MODELER**. **Input:** genome assembly, **Output:** a library containing partially classified consensus sequences of *de-novo* repeat clusters. **Program:** [RepeatModeler](http://www.repeatmasker.org/RepeatModeler.html) and all its many dependencies.

2. **ADD CLASSIFICATIONS FROM THE ONLINE [CENSOR](http://www.girinst.org/censor/) TEXT OUTPUT TO THE REPEAT LIB**. **Input:** A) a library containing partially classified consensus sequences of *de-novo* repeat clusters. B) text file containing the online censor results. **Output:** a library containing more classified consensus sequences of *de-novo* repeat clusters (still some unknow classifications)

3. **SEARCH FOR TEs IN THE GENOME ASSEMBLY WITH REPEAT-MASKER USING THE DENOVO LIB**. **Input:** genome assembly. **Output:** text file with a .out suffix containing classified TE loci with some redundancies. **Program:** [RepeatMasker](http://www.repeatmasker.org/RMDownload.html), which is a dependency of RepeatModeler.

4. **ELIMINATE REDUNDANCIES IN THE REPEAT-MASKER RESULTS**. **Input:** text file with a .out suffix containing classified TE loci with some redundancies, or a directory with several .out files. **Output:** elem_stored.csv files, one per contig. Also generate other files per contig. Puts them in the same folder as the .out files used as input. **Program:** [One Code To Find Them All](http://www.biomedcentral.com/content/pdf/1759-8753-5-13.pdf).

5. **INDEPENDENT SEARCH FOR LTR ELEMENTS BASED ON SECONDARY STRUCTURE**. **Input:** genome assembly, **Output:** text file with loci. **Program:** [LTRharvest](http://www.zbh.uni-hamburg.de/?id=206)

6. **INDEPENDENT SEARCH FOR ELEMENTS BASED ON CODING SEQUENCES**. **Input:** genome assembly, **Output:** text file with loci. **Program:** [TransposonPSI](http://transposonpsi.sourceforge.net/)

7. **ELIMINATE REDUNDANCIES AMONG PROGRAMS**
    1. **Read OneCodeTo... results**. **Input:** Directory containing elem_stored.csv files. **Output:** Dictionary, Integer: num of elements in Dictionary.
    
    2. **Read LTRharvest results and iliminate redundancies chosing the longer match between the programs**. **Input:** Dictionary, Integer: num of elements in Dictionary. **Output:** Dictionary, Integer: num of elements in Dictionary.
    
    3. **Read TransposonPSI results and iliminate redundancies chosing the longer match between the programs**. **Input:** Dictionary, Integer: num of elements in Dictionary. **Output:** Dictionary, Integer: num of elements in Dictionary.
    
8. **PRINT NON-REDUNDANT OUTPUT FILE, ONE PER PROGRAM**. **Input:** Dictionary. **Output:** Three text files.

## Cookbook

### Run RepeatModeler
*RepeatModeler will produce concensus sequeces representing clusters of denovo repeat sequences, partialy classified by RepeatMasker*  

<pre>
from TE import *

make_repeatmodeler_database(name='a_database',
                            input_filename='genome_assembly_file')
# use the BuildDatabase keyword to specify the path to your executable
  
run_repeatmodeler('a_database') 
# use the RepeatModeler keyword to specify the path to your executable
</pre>

The default engine is ncbi and the default lib is eukaryota. RepeatModeler should make a folder in the CWD with the file **'consensi.fa.classified'** in it. It wil write alot of temp files so don't run in Dropbox.

### Run Censor online and add the classifications to your denovo library (optional)

*Copy and paste the results from the webpage into a text file, here named 'Censor_results'.*

<pre>
censor_classifications = parse_online_censor('Censor_results')

print_online_censor(censor_classifications,
                    'censor_classifications_file')

put_censor_classification_in_repeatmodeler_lib('consensi.fa.classified',
                                                censor_classifications,
                                                'consensi.fa.censor')

</pre>

### Run RepeatMasker using the denovo lib
*Some contig names are too long to be parsed in RepeatMasker. However it is possible to replace the names with aliases and have the translations in a file using the first function in this section. It is important to remember to use the aliased genome assembly in the other programs as well, so that redundancies can be resolved.*
  
<pre>
code_sequence_ids('genome_assembly',
                  'names_translations',
                  'coded_genome_assembly',
                  'prefix_for_aliases')

run_repeat_masker('coded_genome_assembly', lib = 'consensi.fa.censor', species=None) 
</pre>
Aliases: if you give 'Mflo' as a prefix, the contig aliases will be 'Mflo_0', 'Mflo_1' ...  
The run_repeat_masker function accepts all the RepeatMasker keywards.
RepeatModeler will write temporary files in a new folder in the CWD so do not run in Dropbox. The output files (most importantly the .out files) will be in the same directory as the input genome assembly.


### Ged rid of redundancies in the .out file
*Two type of rdundancies are possible: 1) within a run, the same locus may have one classification and subsections of it may have other classifications. 2) you may want to make an additional run of RepeatMasker using the eukaryota library instead of the denovo lib. Both types are handled in this stage. You need to put the .out files of all the RepeatMasker runs in one directory and point the function to that directory.*


<pre>
run_OneCodeToFindThemAll('/path/to/RM/.out/files/',
                         'name_of_intermediate_file', 
                         'octfta_output_filename', 
                         'coded_genome_assembly',
                          build_dictionary = 'build_dictionary.pl',
                          octfta = 'one_code_to_find_them_all.pl'
                          )
</pre>
As before, the default path specified in the build_dictionary keyword is the local path on my machine


### Run independent element searches using alternative approaches
*LTRharvest will search for LTR secondary structures and TransposonPSI will do a blastx search against TE protein CDDs*
<pre>
run_LTRharvest('coded_genome_asembly', 
               'name_of_intermideate_file', 
               'ltrharvest_output')                          

run_TransposonPSI('coded_genome_asembly',
                  TPSI = 'perl transposonPSI.pl')
</pre>
The path for TransposonPSI is pointed to by the TPSI keyword.

### Make a non-redundant data structure representing the results of all the searches
*The data structure is a dictionary with the following structure*:
<pre>
TEs = { 'taken': { 'element0': {'ref': {'record': 'the line from the program's output',
                                        'program': 'the program's name'},
                                'contig': 'the contig's name',
                                'start': integer,
                                'end': integer,
                                'length': integer,
                                'lower_tx_level': 'element',
                                'higher_tx_level': 'class or order or family'
                                },

 
                   'element1': {...},
                   
                   
                   'element2': {...},
                   
                   ...


                  },



        'discarded': {...}

      }
</pre>
The internal structure repeats itself within the 'discarded' key. The element number is unique across the taken and discarded elements.

<pre>

\# puting RepeatMasker results as parsed by OneCode... in the data structure
TEs, serial = parse_ocfa_elem_stored('/path/to/RM/.elem_stored.csv'/files/')

\# adding LTRharvest results to the data structure
TEs, serial = integrate_ltrharvest_to_RM_TEs('ltrharvest_output',      
                                             'coded_genome_assembly',  
                                             serial)                      
                                                                       
\# adding TransposonPSI results to the data structure                                                                 
TEs = integrate_TransposonPSI_to_RM_TEs('NameOfInput.TPSI.allHits.chains.bestPerLocus', 
                                        'coded_genome_asembly', 
                                        TEs, 
                                        serial)

</pre>
Redundencies are resolved by taking the longer match across the programs
This data structure is currently the end point of the workflow. I have functions I use to print and plot the data but they are not resonably abstract. I may include them soon. I am likely to include a gff or gtf formater as well.
