# oligopy

# 1. Install oligopy components

Oligopy requires different programs: blast+ (command line blast) and other pip installable packages: primer3, pandas, numpy, joblib, biopython and pyensembl.

The installation of blast+ will require the creation of the desired databases for the comparison of each probe retrieved.

## 1.1. Python environment and dependencies
Make a new envrionment:  
`conda create --name oligopy python=3.8 numpy=1.21 pandas pyarrow joblib`  
`conda activate oligopy`  
Install [Biopython](https://biopython.org/) through pip: `pip install biopython`  
Install [primer3](https://github.com/primer3-org/primer3) through pip: `pip install primer3-py`  
Install [Pyensembl](https://github.com/openvax/pyensembl) through conda: `conda install bioconda::pyensembl`  

## 1.2. Install Blast Command Line Tool
NCBI provides command line standalone BLAST+ programs (based on the NCBI C++ toolkit) as a single compressed package. The package is available for a variety of computer platforms (hardware/operating system combinations) at:  
ftp://ftp.ncbi.nlm.nih.gov/blast/executables/LATEST/  
  
The archives for Linux and Mac OSX are gzip-compressed tar files named using the following convention:  
`ncbi-blast-#.#.#+-CHIP-OS.tar.gz`
Here, the `#.#.#` represents the version number of the current release, CHIP indicates the chipset, and OS indicates the operating system. Equivalent .rpm and .dmg files for Linux and Mac OSX are also available. These archives and their target platforms are listed in the table below.

Installation process from the disk image (.dmg) for Mac OSX and the Red Hat Package Manger (.rpm) for Linux can be found here: http://www.ncbi.nlm.nih.gov/books/NBK279671/ (also explained below).
To install, simply extract the downloaded package after placing it under a desired directory. This can be accomplished by a single tar command, or a combination of gunzip and tar commands.

`tar zxvpf ncbi-blast+#.#.#-x64-linux.tar.gz`

or

`gunzip -d ncbi-blast-#.#.#+-x64-linux.tar.gz`

`tar xvpf ncbi-blast-#.#.#+-x64-linux.tar`

Successful execution of the above commands installs the package and generates a new ncbi-blast-#.#.#+ directory under the working directory selected. This new directory contains the bin and doc subdirectories, as well as a VERSION file. The bin subdirectory contains the programs listed below.
Using the BLAST+ package installed above without configuration could be cumbersome – it requires the installation path to be prefixed to the program call and database specification since the system does not know where to look for the installed program and the specified database. To streamline BLAST searches, two environment variables, PATH and BLASTDB, need to be modified and specified, respectively, to point to the corresponding directories.
Under bash, the following command appends the path to the new BLAST bin directory to the existing PATH setting:

`export PATH=$PATH:$HOME/ncbi-blast-#.#.#+/bin`

The equivalent command under csh is:

`setenv PATH ${PATH}:/home/tao/ncbi-blast-#.#.#+/bin`

The modified $PATH can be examined using echo (added portion underlined):

`echo $PATH`

`/usr/X11R6/bin:/usr/bin:/bin:/usr/local/bin:/opt/local/bin:/home/tao/ncbi-blast-2.2.29+/bin`

To manage available BLAST databases, create a directory to store them:

`mkdir $HOME/blastdb`

Similar approaches described above can be used to set the BLASTDB value under bash:

`export BLASTDB=$HOME/blastdb`

Or under csh to create it anew:

`set BLASTDB=$HOME/blastdb`

A better approach is to have the system automatically set these variables upon login, by modifying the .bash_profile or .cshrc file.
Once they are set, the system knows where to call BLAST programs, and the invoked program will know where to look for the database files. Note that with BLASTDB unspecified, BLAST+ programs only search the working directory, i.e. the directory where BLAST command is issued. For more details about configuring BLAST+, please see http://www.ncbi.nlm.nih.gov/books/NBK279695/.

## 1.3. Make blast database
Currently, oligopy only supports the transcriptome or genome with ENSEMBL format, so that the transcriptome and genome should be downloaded from the ensembl database in fasta format:

Mouse ncRNA and cDNA:  
https://ftp://ftp.ensembl.org/pub/release-87/fasta/mus_musculus/ncrna/  
https://ftp://ftp.ensembl.org/pub/release-87/fasta/mus_musculus/cdna/  

Human ncRNA and cDNA:  
https://ftp://ftp.ensembl.org/pub/release-89/fasta/homo_sapiens/ncrna/  
https://ftp://ftp.ensembl.org/pub/release-89/fasta/homo_sapiens/cdna/  
  
Drosophila ncRNA and cDNA:  
https://ftp.ensembl.org/pub/release-112/fasta/drosophila_melanogaster/ncrna/
https://ftp.ensembl.org/pub/release-112/fasta/drosophila_melanogaster/cdna/

Unzip both files and merge them into a single file that will include all transcriptome in fasta format. This can be done using the bash command:  
`cat Mus_musculus*.fa > Mus_musculus_transcriptome_DB.fa`  
`cat Homo_sapiens*.fa > Homo_sapiens_transcriptome_DB.fa`  
`cat Drosophila_melanogaster*.fa > Drosophila_melanogaster_transcriptome_DB.fa`  

Move the new files to the blastdb folder created in section 1.1. 
In the terminal, move to the blastdb folder and create both databases with the following commands:  
`makeblastdb -in Mus_musculus_transcriptome_DB.fa -parse_seqids -dbtype nucl`  
`makeblastdb -in Homo_sapiens_transcriptome_DB.fa -parse_seqids -dbtype nucl`  
`makeblastdb -in Drosophila_melanogaster_transcriptome_DB.fa -parse_seqids -dbtype nucl`  

## 1.4. Custom sequences

The files HUMAN_NCBI_GENES_retrieved.fasta.masked and MOUSE_NCBI_GENES_retrieved.fasta.masked contain NCBI main isoforms for most of the transcriptome. Oligopy will try to retrieve desired genes from these fasta files (for human or mouse respectively). If not in the file, the sequence can be added manually to these files or oligopy will try to retrieve it automatically from ensembl using pyensembl.

## 1.5. Repeat Masker
For additional instructions see: https://darencard.net/blog/2022-10-13-install-repeat-modeler-masker/  
  
#### Tandem Repeats Finder
Website: https://tandem.bu.edu/trf/home  
From [Github](https://github.com/Benson-Genomics-Lab/TRF/releases/tag/v4.09.1) find the latest release and download the one for your platform.  
Make executable (unix): `chmod -x trf409.macosx`  
Create symlink: `ln -s trf409.macosx trf`

#### Sequence Search Engine
There are a number of options (see RepeatMasker website). We will insall HMMER.  
`brew install hmmer`  
Locate the hmmer folder that contains the nhmmer function and past that path.   

#### RepeatMasker
Download repeat masker from https://www.repeatmasker.org/RepeatMasker/ and install in `/usr/local`.  
`cp RepeatMasker-open-4-#-#.tar.gz /usr/local` If you install it somewhere else, please update the path to the executable in the `oligopy/variables.yaml` file.  
`cd /usr/local`  
`gunzip RepeatMasker-open-4-#-#.tar.gz`  
`tar xvf RepeatMasker-open-4-#-#.tar`  
`cd RepeatMasker`  
`perl ./configure` Choose option 1 to only install the minimal database, downloading will take some time.  
When asked for the Tandem Repeats Finder:   
Enther the made symlink. For me: `/Users/<path to user profile>/Documents/Projects/FISH_general/oligopy/TRF/trf`   
   
When asked for the Sequence Serach Engine:  
Select option 3 for hmmer.  
Locate the hmmer folder that contains the nhmmer function and past that path.  
For me: `/opt/homebrew/Cellar/hmmer/3.4/bin`  
When configured enter `5` for done.  
  
RepeatMasker should now be present in: `/usr/local/RepeatMasker/`  

# 2. Run oligopy: parameters.

Example run with default parameters: `python oligopy.py -query codebookMouse448.xlsx -db Mus_musculus.fa -ncores 12 -db_species mouse -probe_type twist -out Probes`
  
Input can be fasta file with defined sequences, or excel with a column "Gene" and "Tail1" ... "TailN" columns containing the readouts to assign to each gene.
  
`oligopy.py`   
`-h`, show this help message and exit  
`-query`, Input sequences  
`-t`, Minimum probe Tm, default 65    
`-T`, Maximum probe Tm, default 90  
`-size`, Probe length, default 30  
`-ncores`, Number of Cores, default 1  
`-overlap`, Distance between probes, default 2  
`-m`, Minimum probe length, default 26  
`-M`, Maximum probe length, default 32   
`-salt`, Salt concentration for the Tm calculation, default 300  
`-db`, Database: BLASTdatabase mouse or human  
`-start`, Start site to retrieve probes, default 0  
`-end`, End site to retrive probes, default None  
`-mGC`, Minimum GC content, default 0.4  
`-MGC`, Maximum GC content, default 0.6  
`-blast`, Blast Maximum Identity Allowed, default 60  
`-max_probes`, Retrieve maximum of max_probes, default 28  
`-mask`, Start site to retrieve probes. F=False, T=True, default "F"  
`-PNAS`, PNAS, default ["1","2","3","4","5"]  
`-out`, Output file name  
`-Noff`, Allowed number of off target probes with the same off-target match in the same probeset, default 7  
`-db_species`, Choose: human, mouse or drosophila, default None  
`-probe_type`, twist or opool, default "twist"  
