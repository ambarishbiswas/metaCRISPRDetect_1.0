#### metaCRISPRDetect version 1.0:
An ultrafast and highly accurate CRISPR array prediction tool for metagenomes. 

The metaCRISPRDetect program offers rapid identification of CRISPR arrays in short reads with already provided reference CRISPR repeat databases (generated from all publicly available 
archaea and bacterial genomic sequences, as well as CRISPR arrays predicted from several thousands metagenomic assemblies).
 


#### metaCRISPRDetect dependencies:
Please make sure that CRISPRDetect_3.0 (https://github.com/ambarishbiswas/CRISPRDetect_3.0) is properly installed and the program "CRISPRDetect3" is available in the PATH, which should 
satisfy most 3rd party dependencies.

Beside that you only need to install Megahit (version 1.2.9 or higher) or SPAdes (version 3.11.1 or higher) assembler and BBmap suite of tools also known as BBTools.

	Megahit   [Download from https://github.com/voutcn/megahit]
	SPAdes    [Download from https://github.com/ablab/spades]
	BBMap     [Download from https://sourceforge.net/projects/bbmap/]


#### metaCRISPRDetect Syntax:
     
     metaCRISPRDetect -sid ABCD -f test.fasta -o test 
     
#### metaCRISPRDetect examples:
     
     metaCRISPRDetect -sid ABCD -i test.fa -o test 
     metaCRISPRDetect -sid ABCD -s test.fq.gz -o test -read_error_correction 1 -use_ref_repeats_from_assembly 1
     metaCRISPRDetect -sid ABCD -1 forward_fastq_file -2 reverse_fastq_file -o test 


#### metaCRISPRDetect commandline parmeters:

  Options for input sequence(s) in FASTA format :

 	-i/-f	input_fasta_sequence_file	[either gzip compressed (i.e. with extension .gz) or uncompressed FASTA sequence file. Supported extensions are .fa, .fna, .fasta ]
	
	-o/-out_dir	a_folder_name	[A folder with the provided name will be created in the current directory; No "/" or "\" is accepted]
 	
	-sid/-sample_id	TEXT	[This is a compulsory option; typical examples are: NCBI accessions, SRA accessions or any alphanumeric text ]


  Options for input sequence(s) in FASTQ format ONLY :
  
 	-s/-r	input_fastq_sequence_file	[A FASTQ file. Supported extensions are .fq, .fastq with/without gzip compression (i.e. with extension .gz) ]
	-1	forward_fastq_file	[Specify the forward reads FASTQ file of paired-end libraries. Supported extensions are .fq, .fastq with/without gzip compression (i.e. with extension .gz)]
	-2	reverse_fastq_file	[Specify the reverse reads FASTQ file of paired-end libraries. Supported extensions are .fq, .fastq with/without gzip compression (i.e. with extension .gz)]
	--12	interleaved_fastq_file	[Specify a FASTQ file with forward and reverse reads interleaved. Supported extensions are .fq, .fastq with/without gzip compression (i.e. with extension .gz)]


  Options specific for reads in the inputted FASTA/FASQ file :
  
 	-read_error_correction	0/1	[Default is set to 1 (recomended). To skip error correction use 0]
 	-min_read_length	20	[Default value set to 20, Any positive integer is supported]


  Options of assembly of reads in the inputted sequence file :
  
 	-use_ref_repeats_from_assembly	0/1	[Default is set to 0, which means no assembly will be done and CRISPR arrays will only be predicted from the inputted FASTA/FASTQ reads/sequences]
 	-assembler	megahit/spades	[Default set to Megahit; to use SPAdes assembler specify "-assembler spades"]
 	-assembled_contigs	FASTA_seq_file	[specify an existing (multi)FASTA sequence file, which will be used for generating reference CRISPR repeats library. File with extension .gz is supported]


  Options for filtering CRISPR arrays :
  
 	-array_quality_score_cutoff	numeric_value	[default is set to 3, any positive or negative real numeric value is supported]
 	-minimum_no_of_repeats	integer_value	[Default is set to 2, any positive integer >1 is supported]


  Options for additional file/format :
  
 	-create_repeats_fasta_file	0/1	[Default is set to 0, Use "1" to create a multiFASTA sequence file containing the repeats]
 	-create_spacers_fasta_file	0/1	[Default is set to 0, Use "1" to create a multiFASTA sequence file containing the spacers]
 	-create_gff_file	0/1	[Default is set to 0, Use "1" to create a GFF file containing the CRISPR element details]


  Other options :
  
 	-T/-threads	number_of_threads	[By default the program uses 4 threads; Any positive integer is accepted]
 	-q/-quiet	0/1	[Default is set to 0, which shows program step-by-step logs; Use "1" to turn of the logging; Note, a file log.log will still be created in the output_folder]
 	-h/--h/-help/--help		[Shows this help]
 	-v/-version		[Shows version]



If you use metaCRISPRDetect, then please cite:
-----------------------------------------
	Biswas, A., Staals, R. H., Morales, S. E., Fineran, P. C. & Brown, C. M. CRISPRDetect: a flexible algorithm to define CRISPR arrays. BMC Genomics 17, 356 (2016).

For version updates and bug fixes refer to https://github.com/ambarishbiswas/metaCRISPRDetect_1.0 or email ambarishbiswas[at]gmail[dot]com  

 
