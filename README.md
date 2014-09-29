
#README

# TE discovery in a genome assembly

## Overview

These are a set of wrappers compossing a pipline which finds transposable
elements in a genome assembly. The pipline includes the following steps:

1. **MAKE A DENOVO LIB WITH REPEAT-MODELER**. **Input:** genome assembly,
**Output:** a library containing partially classified consensus sequences of
*de-novo* repeat clusters. **Program:**
[RepeatModeler](http://www.repeatmasker.org/RepeatModeler.html) and all its many
dependencies.

2. **ADD CLASSIFICATIONS FROM THE ONLINE
[CENSOR](http://www.girinst.org/censor/) TEXT OUTPUT TO THE REPEAT LIB**.
**Input:** A) a library containing partially classified consensus sequences of
*de-novo* repeat clusters. B) text file containing the online censor results.
**Output:** a library containing more classified consensus sequences of *de-
novo* repeat clusters (still some unknow classifications)

3. **SEARCH FOR TEs IN THE GENOME ASSEMBLY WITH REPEAT-MASKER USING THE DENOVO
LIB**. **Input:** genome assembly. **Output:** text file with a .out suffix
containing classified TE loci with some redundancies. **Program:**
[RepeatMasker](http://www.repeatmasker.org/RMDownload.html), which is a
dependency of RepeatModeler.

4. **ELIMINATE REDUNDANCIES IN THE REPEAT-MASKER RESULTS**. **Input:** text file
with a .out suffix containing classified TE loci with some redundancies, or a
directory with several .out files. **Output:** elem_stored.csv files, one per
contig. Also generate other files per contig. Puts them in the same folder as
the .out files used as input. **Program:** [One Code To Find Them
All](http://www.biomedcentral.com/content/pdf/1759-8753-5-13.pdf).

5. **INDEPENDENT SEARCH FOR LTR ELEMENTS BASED ON SECONDARY STRUCTURE**.
**Input:** genome assembly, **Output:** text file with loci. **Program:**
[LTRharvest](http://www.zbh.uni-hamburg.de/?id=206)

6. **INDEPENDENT SEARCH FOR ELEMENTS BASED ON CODING SEQUENCES**. **Input:**
genome assembly, **Output:** text file with loci. **Program:**
[TransposonPSI](http://transposonpsi.sourceforge.net/)

7. **ELIMINATE REDUNDANCIES AMONG PROGRAMS**
    1. **Read OneCodeTo... results**. **Input:** Directory containing
elem_stored.csv files. **Output:** Dictionary, Integer: num of elements in
Dictionary.

    2. **Read LTRharvest results and eliminate redundancies chosing the longer
match between the programs**. **Input:** Dictionary, Integer: num of elements in
Dictionary. **Output:** Dictionary, Integer: num of elements in Dictionary.

    3. **Read TransposonPSI results and eliminate redundancies chosing the
longer match between the programs**. **Input:** Dictionary, Integer: num of
elements in Dictionary. **Output:** Dictionary, Integer: num of elements in
Dictionary.

8. **PRINT NON-REDUNDANT OUTPUT FILE, ONE PER PROGRAM**. **Input:** Dictionary.
**Output:** Three text files.

## Cookbook

### Run RepeatModeler
*RepeatModeler will produce concensus sequeces representing clusters of denovo
repeat sequences, partialy classified by RepeatMasker*

<pre>
from TE import *

make_repeatmodeler_database(name='a_database',
                            input_filename='genome_assembly_file')
# use the BuildDatabase keyword to specify the path to your executable

run_repeatmodeler('a_database')
# use the RepeatModeler keyword to specify the path to your executable
</pre>

The default engine is ncbi and the default lib is eukaryota. RepeatModeler
should make a folder in the CWD with the file **'consensi.fa.classified'** in
it. It wil write alot of temp files so don't run in Dropbox.

### Run Censor online and add the classifications to your denovo library
(optional)

*Copy and paste the results from the webpage into a text file, here named
'Censor_results'.*

<pre>
censor_classifications = parse_online_censor('Censor_results')

print_online_censor(censor_classifications,
                    'censor_classifications_file')

put_censor_classification_in_repeatmodeler_lib('consensi.fa.classified',
                                                censor_classifications,
                                                'consensi.fa.censor')

</pre>

### Run RepeatMasker using the denovo lib
*Some contig names are too long to be parsed in RepeatMasker. However it is
possible to replace the names with aliases and have the translations in a file
using the first function in this section. It is important to remember to use the
aliased genome assembly in the other programs as well, so that redundancies can
be resolved.*

<pre>
code_sequence_ids('genome_assembly',
                  'names_translations',
                  'coded_genome_assembly',
                  'prefix_for_aliases')

run_repeat_masker('coded_genome_assembly', lib = 'consensi.fa.censor',
species=None)
</pre>
Aliases: if you give 'Mflo' as a prefix, the contig aliases will be 'Mflo_0',
'Mflo_1' ...
The run_repeat_masker function accepts all the RepeatMasker keywards.
RepeatModeler will write temporary files in a new folder in the CWD so do not
run in Dropbox. The output files (most importantly the .out files) will be in
the same directory as the input genome assembly.


### Ged rid of redundancies in the .out file
*Two type of rdundancies are possible: 1) within a run, the same locus may have
one classification and subsections of it may have other classifications. 2) you
may want to make an additional run of RepeatMasker using the eukaryota library
instead of the denovo lib. Both types are handled in this stage. You need to put
the .out files of all the RepeatMasker runs in one directory and point the
function to that directory.*


<pre>
run_OneCodeToFindThemAll('/path/to/RM/.out/files/',
                         'name_of_intermediate_file',
                         'octfta_output_filename',
                         'coded_genome_assembly',
                          build_dictionary = 'build_dictionary.pl',
                          octfta = 'one_code_to_find_them_all.pl'
                          )
</pre>
As before, the default path specified in the build_dictionary keyword is the
local path on my machine


### Run independent element searches using alternative approaches
*LTRharvest will search for LTR secondary structures and TransposonPSI will do a
blastx search against TE protein CDDs*
<pre>
run_LTRharvest('coded_genome_asembly',
               'name_of_intermideate_file',
               'ltrharvest_output')

run_TransposonPSI('coded_genome_asembly',
                  TPSI = 'perl transposonPSI.pl')
</pre>
The path for TransposonPSI is pointed to by the TPSI keyword.

### Make a non-redundant data structure representing the results of all the
searches
*The data structure is a dictionary with the following structure*:
<pre>
TEs = { 'taken': { 'element0': {'ref': {'record': 'the line from the program's
output',
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
The internal structure repeats itself within the 'discarded' key. The element
number is unique across the taken and discarded elements.

<pre>

\# puting RepeatMasker results as parsed by OneCode... in the data structure
TEs, serial = parse_ocfa_elem_stored('/path/to/RM/.elem_stored.csv'/files/')

\# adding LTRharvest results to the data structure
TEs, serial = integrate_ltrharvest_to_RM_TEs('ltrharvest_output',
                                             'coded_genome_assembly',
                                             serial)

\# adding TransposonPSI results to the data structure
TEs = integrate_TransposonPSI_to_RM_TEs('NameOfInput.TPSI.allHits.chains.bestPer
Locus',
                                        'coded_genome_asembly',
                                        TEs,
                                        serial)

</pre>
Redundencies are resolved by taking the longer match across the programs
This data structure is currently the end point of the workflow. I have functions
I use to print and plot the data but they are not resonably abstract. I may
include them soon. I am likely to include a gff or gtf formater as well.


    from Bio.Blast.Applications import NcbitblastnCommandline, NcbipsiblastCommandline
    from Bio.Blast import NCBIXML
    from Bio import SeqIO
    from Bio.Seq import Seq
    from Bio.SeqRecord import SeqRecord
    from Bio import SeqFeature
    from Bio.SeqFeature import SeqFeature, FeatureLocation
    from Bio.Alphabet import IUPAC
    from Bio.Blast import NCBIXML
    import re, os, inspect, subprocess, warnings, sys

ReapeatMasker command line function
===================================


    def code_sequence_ids(in_fasta_file_name, codes_log_filename, out_fasta_file_name, genome_code):
              from Bio import SeqIO
              infile = SeqIO.parse(in_fasta_file_name, 'fasta')
              codes = open(codes_log_filename, 'wt')
              contig_ids_coded = {}
              coded_contigs = []
              count = 1
              for record in infile:
                  contig_ids_coded[genome_code+'_'+str(count)] = record.id
                  record.id = genome_code+'_'+str(count)
                  record.description = ''
                  count += 1
                  coded_contigs.append(record)
              for code in contig_ids_coded.keys():
                  codes.write(code + '\\t' + contig_ids_coded[code] + '\\n')
              SeqIO.write(coded_contigs, out_fasta_file_name, 'fasta')
              codes.close()
              return contig_ids_coded
    
    def run_repeat_masker(query, RepeatMasker = '/home/amir/homeWork/RM_package/RepeatMasker/RepeatMasker', engine=False, parallel=2,
                 slow=False, quick=False, rush=False, nolow=False, noint=False, norna=True, alu=False, div=False,
                 lib=False, cutoff=255, species='eukaryota', is_only=False, is_clip=False, no_is=False, gc=False,
                 gccalc=False, frag=60000, nocut=False, alignments=True, inv=False, lcambig=True,
                 small=False, xsmall=False, poly=False, source=False, html=False, ace=False, gff=False, u=False,
                 xm=False, no_id=True, excln=True, noisy=False):
        
        
        
        frame = inspect.currentframe()
        args, _, _, values = inspect.getargvalues(frame)
        del values['frame']
                        
        # compose the command line 
        
        args_wo_values = ['slow','quick','rush','nolow','noint','norna','alu','is_only','is_clip','no_is'
                          'gccalc','nocut','nocut','alignments','inv','lcambig','small','xsmall','poly',
                          'source','html','ace','gff','u','xm','no_id','excln','noisy']
        
        cline = RepeatMasker + ' '
        for arg in values.keys():
            if not arg in ['RepeatMasker','query']:
                if not arg in args_wo_values and not values[arg] == False:
                    cline = cline + '-' + arg + ' ' + str(values[arg]) + ' '
                elif values[arg] == True:
                    cline = cline + '-' + arg + ' '
    
        cline = cline + query
        
        # execute the command
        print cline
        return os.system(cline), query + '.out' 

# RepeatModeler functions


    def make_repeatmodeler_database(name, input_filename, BuildDatabase='/home/amir/homeWork/RM_package/RepeatModeler/BuildDatabase',
            dir=False, engine='ncbi', batch=False):
        
        frame = inspect.currentframe()
        args, _, _, values = inspect.getargvalues(frame)
        del values['frame']
           
        # compose the command line        
        cline = BuildDatabase + ' '
        for arg in values.keys():
            if not arg in ['BuildDatabase']:
                if not values[arg] == False and not values[arg] == True and not arg == 'input_filename':
                    cline = cline + '-' + arg + ' ' + str(values[arg]) + ' '
        cline = cline + input_filename
        
        # execute the command
        print cline
        return os.system(cline)
    
    def run_repeatmodeler(database, RepeatModeler='perl /home/amir/homeWork/RM_package/RepeatModeler/RepeatModeler', engine='ncbi',
                         species='eukaryota'):
        frame = inspect.currentframe()
        args, _, _, values = inspect.getargvalues(frame)
        del values['frame']
        
        cline = (RepeatModeler + ' -engine ' + engine + ' -database ' + database +
                 ' > ' + database + '.out')
        print cline
        return os.system(cline)
    
    def find_result_dirs_per_genome_code():
        
        """CWD should be a direcrtory containing a set of subdirs with
        RepeatModeler results. This returns a dict with the genome codes
        as keys and the dir name as values"""
        
        RM_results_subdirectories = {}
        for sub in os.walk('.').next()[1]:
            if sub[0:2] == 'RM':
                folder_name = sub
                Genome_code = open('./'+sub+'/round-1/sampleDB-1.fa','r').readlines()[0].split()[1].split('-')[0]
                RM_results_subdirectories[Genome_code]=folder_name
        return RM_results_subdirectories
    #print find_results_dirs_per_genome_code()
    
    # {'Gros': 'RM_6449.WedMay210843002014_Gros'}

#Parse online CENSOR results


    def parse_online_censor(filename):
        censor_classifications = {}
        keep_parsing = True
        lines = open(filename,'r').readlines()
        for line in lines:
            if "Masked Sequence" in line:
                keep_parsing = False
            parse_this_line = True
            NOT = ('[GIRI]', 'Home ', 'Map of Hits', 'SVG', 'Name ')
            for n in NOT:
                if n in line or line[0]=='\t' or line[0]==' ' or line[0]=='\n':
                    parse_this_line = False
            if keep_parsing and parse_this_line:
                l = line.rstrip()
                name = l.split('\t')[0][:-1]
                Class = l.split('\t')[6][:-1]
                score = l.split('\t')[-1]
                if name in censor_classifications.keys():
                    if score > censor_classifications[name]['score']:
                        censor_classifications[name]['score'] = score
                        censor_classifications[name]['Class'] = Class
                else:
                    censor_classifications[name] = {'score': score,
                                                    'Class': Class}
        return censor_classifications
                
    def print_online_censor(censor_classifications, filename):
        import csv
        with open(filename, 'wb') as csvfile:
            linewriter = csv.writer(csvfile, delimiter='\t',
                                    quotechar='|',
                                    quoting=csv.QUOTE_MINIMAL) 
            linewriter.writerow(['Name','calss','score'])
            for name in censor_classifications.keys():
                line = [name, censor_classifications[name]['Class'], censor_classifications[name]['score']]
                linewriter.writerow(line)
                
    def put_censor_classification_in_repeatmodeler_lib(input_filename, censor_classifications, output_filename):
        from Bio import SeqIO
        RM_lib = SeqIO.parse(input_filename, 'fasta')
        RM_CENCOR_lib = []
        for record in RM_lib:
            if 'Unknown' in record.id and record.id in censor_classifications.keys():
                classification = censor_classifications[record.id]['Class']
                record.id = record.id.split('#')[0]+'#'+classification
                record.description = (' ').join(record.description.split(' ')[1:])
                RM_CENCOR_lib.append(record)
            else:
                RM_CENCOR_lib.append(record)
        SeqIO.write(RM_CENCOR_lib,output_filename,'fasta')

# One Code To Find Them All


    def run_OneCodeToFindThemAll(pooled_RM_outputs_dir, #The directory containing repeatmasker .out files or a path to a specific .out file
                                ltr_dict_filename, # name of intermediate file
                                output_filename, 
                                genome_assembly,
                                build_dictionary = '/home/amir/homeWork/RM_package/OCtFtA/build_dictionary.pl',
                                octfta = '/home/amir/homeWork/RM_package/OCtFtA/one_code_to_find_them_all.pl',
                                ):
        cline = build_dictionary+' --rm '+pooled_RM_outputs_dir+' --unknown > '+ ltr_dict_filename
        os.system(cline)
        
        cline = octfta+' --rm '+pooled_RM_outputs_dir+' --ltr '+ltr_dict_filename+' --unknown --fasta '+genome_assembly+' > '+output_filename
        os.system(cline)

# LTRharvest


    def run_LTRharvest(input_filename, index_name, output_name):
        cline = 'gt suffixerator -db '+input_filename+' -indexname '+index_name+' -tis -suf -lcp -des -ssp -sds -dna'
        os.system(cline)
        cline = 'gt ltrharvest -index '+index_name+' -mintsd 5 -maxtsd 100 > '+output_name
        os.system(cline)

# TransposonPSI


    def run_TransposonPSI(input_filename,
                      TPSI = 'perl /home/amir/homeWork/RM_package/TransposonPSI_08222010/transposonPSI.pl'):
        cline = (TPSI+' '+input_filename+' nuc')
        os.system(cline)

# Unite OCTFTA with LTRharvest and TransposonPSI

## Parse OCTFTA elem_stored.csv


    def parse_ocfa_elem_stored(pooled_RM_outputs_dir):
    
        # Get all the '.elem_stored.csv' file names
        from glob import glob
        filenames = glob(pooled_RM_outputs_dir + '*.elem_sorted.csv')
    
        # An empty dict to hold final TE list
        TEs = {'taken': {}, 'discarded': {}}
        serial = 1
    
        # Get all the elements in the OCTFTA output
        for filename in filenames:
            for line in open(filename, 'r').readlines():
                if line[:3] == '###':
                    reference = {'program': 'RMOCFA',
                                 'record': line}
                    contig = line.split('\t')[4]
                    start = line.split('\t')[5]
                    end = line.split('\t')[6]
                    length = line.split('\t')[7]
                    element = line.split('\t')[9]
                    family = line.split('\t')[10]
                    TEs['taken']['element' + str(serial)] = {'ref': reference,
                                                 'contig': contig,
                                                 'start': int(start),
                                                 'end': int(end),
                                                 'length': int(length),
                                                 'lower_tx_level': element,
                                                 'higher_tx_level': family}
                    serial += 1
        return TEs, serial

## Get loci from the LTRharvest output only if they are longer than ones found
with repeatmasker for the same locus


    def integrate_ltrharvest_to_RM_TEs(LTRharvest_output_filename,genome_path, TEs_from_RMOCFA, serial):
        import re
        from Bio import SeqIO
        
        lines = open(LTRharvest_output_filename, 'r').readlines()
        
        contig_names = [record.id for record in SeqIO.parse(genome_path, 'fasta')]
        
        line_count = 0
        for line in lines:
            if not line[0] == '#' and len(line) > 1:
                # correct contig name:
                ## get the true contig name based on the sequence number in the LTRarvest output:
                from Bio import SeqIO
                corrected_sequence_name = None
                line_serial = int(line.rstrip().split('  ')[-1])
                
                try: 
                    corrected_sequence_name = contig_names[line_serial]
                except:
                    raise RuntimeError('Could not find contig for seq number ' + str(line_serial))
                    
                ## Parse the LTRharvest results line
                reference = {'program': 'LTRharvest',
                             'record': line}
                contig = corrected_sequence_name
                start = int(line.split('  ')[0])
                end = int(line.split('  ')[1])
                length = int(line.split('  ')[2])
                lower_tx_level = '?'
                higher_tx_level = 'LTR'
                TE = {'ref': reference,
                      'contig': contig,
                      'start': int(start),
                      'end': int(end),
                      'length': int(length),
                      'lower_tx_level': lower_tx_level,
                      'higher_tx_level': higher_tx_level}
                
                ## Check if the locus is already covered by the repeatmasker results
                ## If it is, check if the ltr hit is longer (then place in taken, and move the rm hit to discarded)
                ## or shorter (then place the ltr hit in discraded)
                placed = False
                for key in TEs_from_RMOCFA['taken'].keys():
                    if (TEs_from_RMOCFA['taken'][key]['start']< start <TEs_from_RMOCFA['taken'][key]['end'] or 
                        TEs_from_RMOCFA['taken'][key]['start']< end <TEs_from_RMOCFA['taken'][key]['end']):
                        ### since it is, keep the longer output (either repeatmasker or LTRharvest)
                        ### use the repeatmasker classification either way
                        ### put the looser in the 'discarded' dictionary
                        if TEs_from_RMOCFA['taken'][key]['length'] < length:
                            #TE['element'] = TEs_from_RMOCFA['taken'][key]['lower_tx_level']
                            #TE['family'] = TEs_from_RMOCFA['taken'][key]['higher_tx_level']
                            TEs_from_RMOCFA['discarded'][key] = TEs_from_RMOCFA['taken'].pop(key, None)
                            TEs_from_RMOCFA['taken']['element'+str(serial)] = TE
                        else:
                            TEs_from_RMOCFA['discarded']['element'+str(serial)] = TE
                        placed = True
                        break
                if not placed:
                    ### Since it is not, add the LTRharvest TE to the 'taken' dict:
                    TEs_from_RMOCFA['taken']['element'+str(serial)] = TE
                serial +=1
                if line_count%100 == 0:
                    print str(line_count)
                line_count += 1
        return TEs_from_RMOCFA, serial

## Get loci from the TransposonPSI output only if they are longer than ones
found with repeatmasker for the same locus


    def integrate_TransposonPSI_to_RM_TEs(TransposonPSI_output_filename,genome_path, TEs_from_RMOCFA, serial):
        import re
        
        lines = open(TransposonPSI_output_filename, 'r').readlines()
        
        line_count = 0
        for line in lines:
            if line[0] == '#':
                    
                ## Parse the TransposonPSI results line
                reference = {'program': 'TransposonPSI',
                             'record': line}
                contig = line.split('\t')[3]
                start = int(line.split('\t')[4].split('-')[0])
                end = int(line.split('\t')[4].split('-')[1])
                length = end-start+1
                lower_tx_level = '?'
                higher_tx_level = line.split('\t')[1]
                TE = {'ref': reference,
                      'contig': contig,
                      'start': int(start),
                      'end': int(end),
                      'length': int(length),
                      'lower_tx_level': lower_tx_level,
                      'higher_tx_level': higher_tx_level}
                
                ## Check if the locus is already covered by previous results
                placed = False
                for key in TEs_from_RMOCFA['taken'].keys():
                    if (TEs_from_RMOCFA['taken'][key]['start']< start <TEs_from_RMOCFA['taken'][key]['end'] or 
                        TEs_from_RMOCFA['taken'][key]['start']< end <TEs_from_RMOCFA['taken'][key]['end']):
                        ### since it is, keep the longer output 
                        ### put the looser in the 'discarded' dictionary
                        if TEs_from_RMOCFA['taken'][key]['length'] < length:
                            TEs_from_RMOCFA['discarded'][key] = TEs_from_RMOCFA['taken'].pop(key, None)
                            TEs_from_RMOCFA['taken']['element'+str(serial)] = TE
                        else:
                            TEs_from_RMOCFA['discarded']['element'+str(serial)] = TE
                        placed = True
                        break
                if not placed:
                    ### Since it is not, add the TransposonPSI TE to the 'taken' dict:
                    TEs_from_RMOCFA['taken']['element'+str(serial)] = TE
                serial +=1
                if line_count%100 == 0:
                    print str(line_count)
                line_count += 1
        return TEs_from_RMOCFA


    
