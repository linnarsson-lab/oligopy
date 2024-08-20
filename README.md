# oligopy
Oligopy is a parallelized pure python program to design smFISH probes.  
Oligopy can add tails that encode a binary barcode, as well as PCR primers for probe amplification. The probes either have a U-shape if they are PCR amplified. Or an L-shape if they are synthesized chemically, where the RNA binding part is located at the 5' end. After probes are designed Excel order forms will be made for IDT or Twist Bioscience, so that probes can be directly ordered. 

# 1. Install oligopy components

Oligopy requires different programs: blast+ (command line blast) and other pip installable packages: primer3, pandas, numpy, joblib, biopython and pyensembl.  
The installation of blast+ will require the creation of the desired databases for the comparison of each probe retrieved.

## 1.1. Python environment and dependencies
Make a new environment:  
`conda create --name oligopy python=3.9 numpy pandas joblib`  
`conda activate oligopy`  
Install [Biopython](https://biopython.org/) through pip: `pip install biopython`  
Install [primer3](https://github.com/primer3-org/primer3) through pip: `pip install primer3-py`  
Install [Pyensembl](https://github.com/openvax/pyensembl) through pip: `pip install pyensembl`  
Install openpyxl through pip: `pip install openpyxl`  
Or all together: `pip install biopython primer3-py pyensembl openpyxl`  
Install pytables through conda: `conda install -c conda-forge pytables`  
(Do not use a higher Python version and do use Pip where indicated.)

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

A better approach is to have the system automatically set these variables upon login, by modifying the .bash_profile, .bashrc or .cshrc file.  
Example with bashrc:  
`nano ~/.bashrc`  
Add the path to the Blast bin to the end of the file:  
`export PATH=$PATH:<path to ncbi folder>/ncbi-blast-2.15.0+/bin`  
Close and source the bashrc again: `source ~/.bashrc`  
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

## 1.4. Gene sequences
Oligopy has 3 ways to find the transcript reference sequence for the requested gene to design the probes against. First it will look for the fasta file given by the user through the `-reference_sequence` argument. If this is not given or not found, it will use the given species name through the `-db_species` argument, to find pre-made reference sequences that are stored in the `Oligopy/reference_transcripts` folder. If this file is also not found, Oligopy will use pyensemble to fetch the longest transcript present in Ensembl. Oligopy will also do this if one or more genes could not be found in the file.  
Currently there are reference sequences made for:
    - homo sapiens: Made by Alejandro, repeats are masked with RepeatMasker.  
    - mus musculus: Made by Alejandro, repeats are masked with RepeatMasker.  
    - drosophila melanogaster: Ensembl canonical transcripts.  
If you add new files to the `Oligopy/reference_transcripts` folder, preferably use the latin species name with a lower dash, like: homo_sapiens.  
Make sure you only have one reference file per species. If you want to swap between multiple references, you can provide the file directly via the `-reference_sequence` argument.  
For more information on the ensemble canonical transcript variant see: https://www.ensembl.org/info/genome/genebuild/canonical.html  

## 1.5. Pyensembl
You need to install the correct pyensembl database for the species you want to use. 
In a terminal with the `oligopy` environment activated run: `pyensembl install --release 110 --species drosophila_melanogaster`  
Where you replace the release number and species for the desired release and species.  
For supported species and suggested releases see: `variables.ini`.  
Follow the same instructions when you get an message like: `ValueError: GTF database needs to be created, run: pyensembl install --release 110 --species drosophila_melanogaster`  
   
Development note:  
To add new species you need to adapt the source code and install the new species through pyensembl.  
First check which species and which releases are available [here](https://github.com/openvax/pyensembl/issues/296)  
Then add them in the conda environment by running: `pyensembl install --release 103 110 --species drosophila_melanogaster`  
Add the species name and release to `variables.ini`  

## 1.6. Repeat Masker
Only needed when you want to mask repetitive sequences in your own reference transcripts or if you want to use the `-mask` option.  
Not a requirement for use with default Mouse and Human transcriptome, because the masked fasta files are already made.
  
Install from: https://www.repeatmasker.org/ and paste executable path in the variables.ini file.
For additional instructions see: https://darencard.net/blog/2022-10-13-install-repeat-modeler-masker/  
Nevertheless, installation can be challenging, as an alternative you can use RepeatMasker from this docker container https://github.com/Dfam-consortium/TETools  
However, it will not be possible to use the `-mask` option in that case.  

#### Tandem Repeats Finder
Website: https://tandem.bu.edu/trf/home  
From [Github](https://github.com/Benson-Genomics-Lab/TRF/releases/tag/v4.09.1) find the latest release and download the one for your platform.  
Make executable (unix): `chmod -x trf409.macosx`  
Create TRF symlink: `ln -s trf409.macosx trf`

#### Sequence Search Engine
There are a number of options (see RepeatMasker website). We will install HMMER.  
`brew install hmmer`  
Locate the hmmer folder: `where nhmmer`, save this path.   

#### RepeatMasker
Download repeat masker from https://www.repeatmasker.org/RepeatMasker/ and install in `/usr/local`.  
`cp RepeatMasker-open-4-#-#.tar.gz /usr/local` If you install it somewhere else, please update the path to the executable in the `oligopy/variables.yaml` file.  
`cd /usr/local`  
`gunzip RepeatMasker-open-4-#-#.tar.gz`  
`tar xvf RepeatMasker-open-4-#-#.tar`  
`cd RepeatMasker`  
`perl ./configure` Choose option 1 to only install the minimal database, downloading will take some time.  
If you want to use different databases for different species please download them from here: https://github.com/Dfam-consortium/FamDB  
When asked for the Tandem Repeats Finder:   
Enter the TRF symlink. Example: `/Users/<path to user profile>/Documents/Projects/FISH_general/oligopy/TRF/trf`   
   
When asked for the Sequence Search Engine:  
Select option 3 for hmmer.  
Locate the hmmer folder that contains the nhmmer function and past that path.  
For me: `/opt/homebrew/Cellar/hmmer/3.4/bin`  
When configured enter `5` for done.  
  
RepeatMasker should now be present in: `/usr/local/RepeatMasker/`  
From that folder you can run `./RepeatMasker -h`  
Add the executable path to the `variables.ini` file.

# 2. Run oligopy: parameters.

Example run with default parameters: `python oligopy.py -query codebookMouse448.xlsx -db Mus_musculus.fa -ncores 12 -db_species mouse -probe_type twist -out Probes`
  
Input can be <s>a fasta file with defined sequences</s>, or excel with a column "Gene" and "Tail1" ... "TailN" columns containing the readouts to assign to each gene.
  
`oligopy.py`   
`-h`, show this help message and exit  
`-query`, Input sequences  
`-t`, Minimum probe melting temperature, default 65    
`-T`, Maximum probe melting temperature, default 90  
`-size`, Probe length, default 30  
`-ncores`, Number of CPU cores, default 1  
`-overlap`, Distance between probes, default 2  
`-m`, Minimum probe length, default 26  
`-M`, Maximum probe length, default 32   
`-salt`, Sodium salt concentration for the melting temperature calculation, default 330  
`-db`, Database: BLASTdatabase mouse or human  
`-start`, Start site to retrieve probes, default 0  
`-end`, End site to retrive probes, default None  
`-mGC`, Minimum GC content, default 0.4  
`-MGC`, Maximum GC content, default 0.6  
`-blast`, Blast Maximum Identity Allowed, default 60  
`-max_probes`, Maximum number of probes to retrieve per gene, default 30  
`-mask`, Mask repetitive sequences. Requires RepeatMasker. Not functional at the moment. F=False, T=True, default "F"  
`-PNAS`, PNAS rules to apply, default ["1","2","3","4","5"]  
`-out`, Output file name  
`-Noff`, Allowed number of off target probes with the same off-target match in the same probeset, default 7  
`-db_species`, Species name. Should be supported by pyensembl, default None  
`-probe_type`, Sequence order forms in format for "twist", "opool" or "opool_amp". Twist = Fw primer-readout 1:3 - Probe - Readout 4:6 - Rv primer. Opool = Probe - Readout 1:6. Opool_amp = Fw primer - Probe - Readout 1:6 - Rv primer, default "twist"  
`-max_probes_overlapping`, Maximum number of probes allowed to overlap, default 4  
`-min_probes`, Minimum probes to activate high overlapping mode", default 10  
`-assign_tails`, Add tails to the probes when input is .xlxs. T=True, F=False", default "T"  
`-cleanup`, If `T` Delete intermediate files. F for False, default "T"  
`-reference_sequence`, Filename of fasta file with reference transcripts of each genes, default = "None"