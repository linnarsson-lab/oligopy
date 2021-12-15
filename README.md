# oligopy

## 1. Install oligopy components

Oligopy requires different programs: blast+ (command line blast) and other pip installable packages: primer3, pandas, numpy, joblib, biopython and pyensembl.

The installation of blast+ will require the creation of the desired databases for the comparison of each probe retrieved.

# 1.1. Install Blast Command Line Tool
NCBI provides command line standalone BLAST+ programs (based on the NCBI C++ toolkit) as a single compressed package. The package is available for a variety of computer platforms (hardware/operating system combinations) at:
ftp://ftp.ncbi.nlm.nih.gov/blast/executables/LATEST/
The archives for Linux and Mac OSX are gzip-compressed tar files named using the following convention:
ncbi-blast-#.#.#+-CHIP-OS.tar.gz
Here, the #.#.# represents the version number of the current release, CHIP indicates the chipset, and OS indicates the operating system. Equivalent .rpm and .dmg files for Linux and Mac OSX are also available. These archives and their target platforms are listed in the table below.

Installation process from the disk image (.dmg) for Mac OSX and the Red Hat Package Manger (.rpm) for Linux can be found here: http://www.ncbi.nlm.nih.gov/books/NBK279671/ (also explained below).
To install, simply extract the downloaded package after placing it under a desired directory. This can be accomplished by a single tar command, or a combination of gunzip and tar commands.

$ tar zxvpf ncbi-blast+#.#.#-x64-linux.tar.gz
or
$ gunzip -d ncbi-blast-#.#.#+-x64-linux.tar.gz
$ tar xvpf ncbi-blast-#.#.#+-x64-linux.tar

Successful execution of the above commands installs the package and generates a new ncbi-blast-#.#.#+ directory under the working directory selected. This new directory contains the bin and doc subdirectories, as well as a VERSION file. The bin subdirectory contains the programs listed below.
Using the BLAST+ package installed above without configuration could be cumbersome – it requires the installation path to be prefixed to the program call and database specification since the system does not know where to look for the installed program and the specified database. To streamline BLAST searches, two environment variables, PATH and BLASTDB, need to be modified and specified, respectively, to point to the corresponding directories.
Under bash, the following command appends the path to the new BLAST bin directory to the existing PATH setting:

$ export PATH=$PATH:$HOME/ncbi-blast-#.#.#+/bin
The equivalent command under csh is:
$ setenv PATH ${PATH}:/home/tao/ncbi-blast-#.#.#+/bin
The modified $PATH can be examined using echo (added portion underlined):
$ echo $PATH
/usr/X11R6/bin:/usr/bin:/bin:/usr/local/bin:/opt/local/bin:/home/tao/ncbi-blast-2.2.29+/bin
To manage available BLAST databases, create a directory to store them:
$ mkdir $HOME/blastdb
Similar approaches described above can be used to set the BLASTDB value under bash:
$ export BLASTDB=$HOME/blastdb
Or under csh to create it anew:
set BLASTDB=$HOME/blastdb
A better approach is to have the system automatically set these variables upon login, by modifying the .bash_profile or .cshrc file.
Once they are set, the system knows where to call BLAST programs, and the invoked program will know where to look for the database files. Note that with BLASTDB unspecified, BLAST+ programs only search the working directory, i.e. the directory where BLAST command is issued. For more details about configuring BLAST+, please see http://www.ncbi.nlm.nih.gov/books/NBK279695/.
# 1.2. Make blast database
Currently, oligopy only supports the transcriptome or genome with ENSEMBL format, so that the transcriptome and genome should be downloaded from the ensembl database in fasta format:

Mouse ncRNA and cDNA:
ftp://ftp.ensembl.org/pub/release-87/fasta/mus_musculus/ncrna/
ftp://ftp.ensembl.org/pub/release-87/fasta/mus_musculus/cdna/

Human ncRNA and cDNA:
ftp://ftp.ensembl.org/pub/release-89/fasta/homo_sapiens/ncrna/
ftp://ftp.ensembl.org/pub/release-89/fasta/homo_sapiens/cdna/

Unzip both files and merge them into a single file that will include all transcriptome in fasta format. This can be done using the bash command:
cat Mus_musculus*.fa > Mus_musculus_transcriptome_DB.fa
cat Homo_sapiens*.fa > Homo_sapiens_transcriptome_DB.fa

Move the new files to the blastdb folder created in section 1.1. 
In the terminal, move to the blastdb folder and create both databases with the following commands:
makeblastdb -in Mus_musculus_transcriptome_DB.fa -parse_seqids -dbtype nucl
makeblastdb -in Homo_sapiens_transcriptome_DB.fa -parse_seqids -dbtype nucl

## 2. Run oligopy: parameters.

Example run with default parameters: python oligopy.py -query codebookMouse448.xlsx -db Mus_musculus.fa -ncores 12 -db_species mouse -probe_type twist -out Probes

Input can be fasta file with defined sequences, or excel with a column "Gene" and "Tail1" ... "TailN" columns containing the readouts to assign to each gene.

oligopy.py [-h] [-query Input sequences] [-t Minimum probe Tm]
                    [-T Maximum probe Tm] [-size Probe length]
                    [-ncores Number of Cores]
                    [-overlap Distance between probes]
                    [-m Minimum probe length] [-M Maximum probe length]
                    [-salt Salt concentration for the Tm calculation]
                    [-db Database: BLASTdatabase mouse or human]
                    [-start Start site to retrieve probes]
                    [-end End site to retrive probes]
                    [-mGC Minimum GC content] [-MGC Maximum GC content]
                    [-blast Blast Maximum Identity Allowed]
                    [-max_probes Retrieve maximum of max_probes, default 28]
                    [-mask Start site to retrieve probes] [-PNAS PNAS]
                    [-out Output file name]
                    [-Noff Allowed number of off target probes with the same off-target match in the same probeset]
                    [-db_species Choose: human, mouse or rat]
                    [-probe_type twist or opool]
