% [WorkShop de Verão Bioinformática](http://github.com/genomika/summercourse/)
% __Introdução ao Linux__
% _(updated 21-01-2014)_

<!-- COMMON LINKS HERE -->

# Getting around in the shell

## Important keyboard shortcuts

Type as little and as accurately as possible by using keyboard shortcuts!

### Tab key completion

The Tab key is your best friend! Hit the Tab key once or twice - it's almost always magic! Hitting Tab invokes "shell completion", instructing the shell to try to guess what you're doing and finish the typing for you. On most modern Linux shells, Tab completion will:
-   complete file or directory names up to any ambiguous part (single Tab)
-   if nothing shows up, there is no unambiguous match
-   display all possible completions (Tab twice)
-   you then decide where to go next
-   work for shell commands too (like rsync or chmod)


### Up arrow
Use "up arrow" to retrieve any of the last 500 commands you've typed, going backwards through your history. You can then edit them and hit Enter (even in the middle of the command) and the shell will use that command. The down arrow "scrolls" forward from where you are in the command history.


### Ctrl-a, Ctrl-e
You can use Ctrl-a (holding down the "control" key and "a") to jump the cursor right to the beginning of the line. The omega to that alpha is Ctrl-e, which jumps the cursor to the end of the line. Arrow keys work, and Ctrl-arrow will skip by word forward and backward.

### Wildcards and special file names

The shell has shorthand to refer to groups of files by allowing wildcards in file names.
* (asterisk) is the most common filename wildcard. It matches "any length of any characters".

* Other useful ones are brackets ( [ ] ) to allow for any character in the list of characters between the brackets. And you can use a hyphen ( - ) to specify a range of characters

For example:
-   **ls *.bam** – lists all files in the current directory that end in .bam
-   **ls [A-Z]*.bam** – does the same, but only if the first character of the file is a capital letter
-   **ls [ABab]*.bam** – lists all .bam files whose 1st letter is A, B, a or b.

Three special file names:
-   . (single period) means "this directory".
-   .. (two periods) means "directory above current." So ls -l .. means "list contents of the parent directory."
-   ~ (tilde) means "my home directory".


### Environment variables

Environment variables are just like variables in a programming language (in fact bash is a complete programming language), they are "pointers" that reference data assigned to them. In bash, you assign an environment variable as shown below:

    export varname="Some value, here it's a string"


**NOTE:** Be careful – do not put spaces around the equals sign when assigning environment variable values. Also, always use double quotes if your value contains (or might contain) spaces.

You set environment variables using the bare name (varname above).
You then refer to or evaluate an environment variables using a dollar sign ( $ ) before the name:

    echo $varname

The export keyword when you're setting ensures that any sub-processes that are invoked will inherit this value. Without the export only the current shell process will have that variable set.

Use the env command to see all the environment variables you currently have set

## Standard streams

Every command and Linux program has three "built-in" streams: standard input, standard output and standard error.


![Terminal Example](https://wikis.utexas.edu/download/attachments/66698131/std_streams.png?version=1&modificationDate=1400467628000&api=v2)

### redirecting output
-   To take the standard output of a program and save it to a file, you use the ```>``` operator
-   a single ```>``` overwrites any existing target; a double ```>>``` appends to it
-   since standard output is stream #1, this is the same as ```>1```
-   To redirect the standard error of a program you must specify its stream number using```2>```
-   To redirect standard output and standard error to the same place, use the syntax ```2>&1```

To see the difference between standard output and standard error try these commands:

    # redirect a long listing of your $HOME directory to a file
    ls -la $HOME > cmd.out
    # look at the contents -- you'll see just files
    cat cmd.out
    # this command gives an error because the target does not exist
    ls -la bad_directory
    # redirect any errors from ls to a fil
    ls -la bad_directory 2> cmd.out
    # look at the contents -- you'll see an error message
    cat cmd.out
    # now redirect both error and output streams to the same place
    ls -la bad_directory $HOME > cmd.out
    # look at the contents -- you'll see both an error message and files
    cat cmd.out


It is easy to not notice the difference between standard output and standard error when you're in an interactive Terminal session – because both outputs are sent to the Terminal. But they are separate streams, with different meanings. When running batch programs and scripts you will want to manipulate standard output and standard error from programs appropriately.


## A Janela de Terminal
- Macs e Linux possuem programas de terminal nativos - encontre-os em sua máquina
- Usuários Windows vão precisar de ajuda, opções a seguir:
    - [Putty](http://www.chiark.greenend.org.uk/~sgtatham/putty/download.html) "Putty"
    - [Git-bash](http://msysgit.github.io/) "GitBash"
    - [Cygwin](http://www.cygwin.com/) "Cygwin"

## SSH

ssh é um programa executável que roda em seu computador local e permite que você se conecte de forma segura a um computador remoto. Em Macs, Linux e Windows (Git-Bash out Cygwin), você pode rodar de uma tela de terminal. Responda yes a questão de segurança do prompt do SSH.

##### Acessando servidor remoto via SSH.
    ssh seuusuario@172.16.225.107

- Se você está usando Putty como seu terminal no Windows:
    -  Double-click the Putty.exe icon
    -  In the PuTTY Configuration window
    -  make sure the Connection type is SSH
    -  enter stampede.tacc.utexas.edu for Host Name
    -  click Open button
    -  answer Yes to the SSH security question
    -  In the PuTTY terminal
    -  enter your TACC user id after the login as: prompt, then Enter


## Shell Bash

You're now at a command line! It looks as if you're running directly on the remote computer, but really there are two programs communicating: your local Terminal and the remote Shell. There are many shell programs available in Linux, but the default is bash (Bourne-again shell). The Terminal is pretty "dumb" – just sending your typing over its secure sockets layer (SSL) connection to TACC, then displaying the text sent back by the shell. The real work is being done on the remote computer, by programs called by the bash shell.

![Terminal Example](https://wikis.utexas.edu/download/attachments/66696867/terminal.png?version=1&modificationDate=1400007439000&api=v2)

##### Preparando seu ambiente de trabalho.

First create a few directories and links we will use (more on these later).

**NOTE:** You can copy and paste these lines from the code block below into your Terminal window. Just make sure you hit "Enter" after the last line.

    cd
    ln -s -f $SCRATCH scratch
    ln -s -f $WORK work
    ln -s -f /corral-repl/utexas/BioITeam

- $WORK and $SCRATCH are TACC environment variable that refer to your work and scratch file system areas.
- The ln -s command creates a symbolic link, a shortcut the the linked file or directory.
- here the link targets are your work and scratch file system areas
- having these link shortcuts will help when you want to copy files to your work or scratch, and when you navigate the TACC file system using a remote SFTP client
- always change directory (cd) to the directory where we want the links created before executing ln -s 
- here we want the links under your home directory (cd with no arguments)

##### Want to know where a link points to? Use ls with the -l (long listing) option.
    ls -l

Set up a $HOME/local/bin directory and link some scripts there that we will use a lot in the class.

    mkdir -p $HOME/local/bin
    cd $HOME/local/bin
    ln -s -f /corral-repl/utexas/BioITeam/bin/launcher_creator.py
    ln -s -f /work/01063/abattenh/local/bin/cutadapt
    ln -s -f /work/01063/abattenh/local/bin/samstat

-  The mkdir command creates a new directory. The -p option says to create intermediate directories if needed (like local here).
-  here we're creating a $HOME/local/bin directory where we'll put some programs used in the course
-  $HOME is an environment variable set by TACC that refers to your home directory.
-  The ln -s command creates a symbolic link, a shortcut the the linked file or directory.
-  here the link targets are programs we want – instead of copying the programs, we just link to them
-  always change directory (cd) to the directory where we want the links created before executing ln -s
-  here we want the links in $HOME/local/bin

Want to know more about a Linux command? Type the command name then the --help option. For example, with mkdir:

    mkdir --help

This won't work all the time, but it's your best 1st choice for help.

Now execute the lines below to set up a login script, called .profile_user.

Whenever you login via an interactive shell as you did above, a well-known script is executed by the shell to establish your favorite environment settings. We've set up a common profile for you to start with that will help you know where you are in the file system and make it easier to access some of our shared resources. To set up this profile, do the steps below:

    cd
    cp /corral-repl/utexas/BioITeam/core_ngs_tools/common/stampede_dircolors .dircolors
    cp /corral-repl/utexas/BioITeam/core_ngs_tools/common/core_ngs_profile .profile_user
    chmod 600 .profile_user


-  The chmod 600 .profile_user command marks the file as readable and writable only by you. The .profile_user script file will not be executed unless it has these exact permissions settings.
-  The well-known filename is .profile_user (or .profile on some systems), which is specific to the bash shell.

Since .profile_user is executed when you login, to ensure it is set up properly you should first log off stampede like this:

    exit
    
Then log back in to stampede.tacc.utexas.edu. This time your .profile_user will be executed and you should see a new shell prompt:
    
    stamp:~$

The great thing about this prompt is that it always tells you where you are, which avoids having to issue the pwd (present working directory) command all the time. Execute these commands to see how the prompt reflects your current directory. (Don't just copy-and-paste here because we've included the prompt.)

    stamp:~$ mkdir -p tmp/a/b/c
    stamp:~$ cd tmp/a/b/c
    stamp:~/tmp/a/b/c$

The prompt now tells you you are in the c sub-directory of the b sub-directory of the a sub-directory of the tmp sub-directory of your home directory ( ~ ).

**NOTE:**  The tilde character ( ~ ) is a shortcut that means "home directory". We'll see more of it later.

Your profile has also installed nice directory colors, which you see when you list your home directory:

    cd  
    ls

So why don't you see the .profile_user file you copied to your home directory? Because all files starting with a period ("dot files") are hidden by default. To see them add the -a (all) option to ls:

    ls -a

To see even more detail, including file type and permissions and symbolic link targets, add the -l (long listing) switch:

    ls -la

### Details about your login profile

We list its content to the Terminal with the cat (concatenate files) command that simply reads a file and writes each line of content to standard output (here, your Terminal):

    cat .profile_user
    

**NOTE:**   The cat command just echos the entire file's content, line by line, without pausing, so should not be used to display large files. Instead, use a "pager" (like more or less) or look at parts of the file with **head** or **tail**.

You'll see the following (you may need to scroll up a bit to see the beginning):

    #!/bin/basH
    # Change the command line prompt to contain the current directory path
    if [ "$TACC_SYSTEM" == "stampede" ]; then
        PS1='stamp:\w$ '
    else
        PS1='lstar:\w$ '
    fi
    # Try to ensure all created files can be read/writtin by group members
    umask 002
    # Make common, useful software always available
    module load python; module load launcher
    # Set the default project allocation for launcher_creator.py
    export ALLOCATION=genomeAnalysis
    # Environment variables for useful locations
    export BI=/corral-repl/utexas/BioITeam
    export CLASSDIR="$BI/core_ngs_tools"
    # Add current directory and $HOME/local/bin to PATH
    export PATH=.:$HOME/local/bin:$PATH
    # Use yellow for directories, not that horrible blue
    dircolors .dircolors > /dev/null


So what does the common profile file do? Several things. Let's look at a few of them.

#### she-bang

The first line is the "she-bang". It tells the shell what program should execute this file – in this case, bash itself – even though the expression is inside a shell comment (denoted by the # character).

    #!/bin/bash

#### environment variables

The profile also sets an environment variable named BI to point to the shared directory: **/corral-repl/utexas/BioITeam**, and another environment variable named **CLASSDIR** to point to the specific sub-directory for our class.

    # Environment variables for useful locations
    export BI=/corral-repl/utexas/BioITeam
    export CLASSDIR="$BI/core_ngs_tools"

Environment variables are like variables in a programming language like python or perl (in fact bash is a complete programming language). They have a name (like BI above) and a value (the value for BI is the pathname **/corral-repl/utexas/BioITeam**).

#### shell completion

You can use these environment variables to shorten typing, for example, to look at the contents of the shared BioITeam directory as shown below, using the magic Tab key to perform shell completion.

    # hit Tab once after typing $BI/ to expand the environment variable
    ls $BI/
    # now hit Tab twice to see the contents of the directory
    ls /corral-repl/utexas/BioITeam/
    # now type "co" and hit Tab again
    ls /corral-repl/utexas/BioITeam/co
    # your command line should now look like this
    ls /corral-repl/utexas/BioITeam/core_nge_tools/
    # now type "m" and one Tab
    ls /corral-repl/utexas/BioITeam/core_nge_tools/m
    # now just type one Tab
    ls /corral-repl/utexas/BioITeam/core_nge_tools/misc/
    # the shell expands as far as it can unambiguously, so your command line should look like thi
    ls /corral-repl/utexas/BioITeam/core_nge_tools/misc/small
    # type a period (".") then hit Tab twice again -- you're narrowing down the choices
    ls /corral-repl/utexas/BioITeam/core_nge_tools/misc/small.
    # finally, hit Tab twice to see possible completions now -- you should see two filenames

#### Important Tip -- the Tab key is your BFF!
-   The Tab key is one of your best friends in Linux. Hitting it invokes "shell completion", which is as close to magic as it gets!
-   Tab once will expand the current command line contents as far as it can unambiguously. 
-   if nothing shows up, there is no unambiguous match
-   Tab twice will give you a list of everything the shell finds matching the current command line.
-   you then decide where to go next


#### extending the $PATH

When you type a command name the shell has to have some way of finding what program to run. The list of places (directories) where the shell looks is stored in the $PATH environment variable. You can see the entire list of locations by doing this:

    echo $PATH

As you can see, there are a lot of locations on the **$PATH**. That's because when you load modules at TACC (such as the module load lines in the common profile), that mechanism makes the programs available to you by putting their installation directories on your $PATH. We'll learn more about modules shortly.

Here's how the shared profile adds your **$HOME/local/bin** directory to the location list – recall that's where we linked some programs we'll use – along with a special dot character ( . ) that means "here", or "whatever the current directory is".

    # Add current directory and $HOME/local/bin to PATH 
    export PATH=.:$HOME/local/bin:$PATH

#### setting up the friendly command prompt

The complicated looking if statement near the top of your profile is checking whether you're on stampede or lonestar (this .profile_user works on both), and setting up your friendly shell prompt so that it includes the current working directory. This is done by setting the special PS1 environment variable and including a special \w directive that the shell knows means "current directory".

    # Change the command line prompt to contain the current directory path
    if [ "$TACC_SYSTEM" == "stampede" ]; then
        PS1='stamp:\w$ '
    else
        PS1='lstar:\w$ '
    fi


### Download from a link – wget

Well, you don't have a desktop at TACC to "Save as" to, so what to do with a link? The wget program knows how to access web URLs such as http, https and ftp.

#### wget

Get ready to run wget from the directory where you want to put the data. Don't press Enter after the wget command – just put a space after it.

    cd $WORK/archive/original/2014_05.core_ngs
    wget

Here are two web links:
-   [http://web.corral.tacc.utexas.edu/BioITeam/yeast_stuff/Sample_Yeast_L005_R1.cat.fastq.gz] (http://web.corral.tacc.utexas.edu/BioITeam/yeast_stuff/Sample_Yeast_L005_R1.cat.fastq.gz)
-   [http://web.corral.tacc.utexas.edu/BioITeam/yeast_stuff/Sample_Yeast_L005_R2.cat.fastq.gz] (http://web.corral.tacc.utexas.edu/BioITeam/yeast_stuff/Sample_Yeast_L005_R2.cat.fastq.gz)

Right click on the 1st link in your browser, then select "Copy link location" from the menu. Now go back to your Terminal. Put your cursor after the space following the wget command then either right-click, or Paste. The command line to be executed should look like this:

    wget http://web.corral.tacc.utexas.edu/BioITeam/yeast_stuff/Sample_Yeast_L005_R1.cat.fastq.gz

Now press Enter to get the command going. Repeat for the 2nd link. Check that you now see the two files (ls).


#### Copy from a corral location - cp or rsync

Suppose you have a corral allocation where your organization keeps its data, and that the sequencing data has been downloaded there. You can use various Linux commands to copy the data locally from there to your $SCRATCH area.

#### cp

The cp command copies one or more files from a local source to a local destination. It has the most common form:

    cp [options] <source file 1> <source file 2> ... <destination directory>

Make a directory in your scratch area and copy a single file to it. The trailing slash ( / ) on the destination says it is a directory.

    mkdir -p $SCRATCH/data/test1
    cp  /corral-repl/utexas/BioITeam/web/tacc.genomics.modules  $SCRATCH/data/test1/
    ls $SCRATCH/data/test1

Copy a directory to your scratch area. The -r argument says "recursive".

    cds
    cd data
    cp -r /corral-repl/utexas/BioITeam/web/general/ general/


#### local rsync

The rsync command is typically used to copy whole directories. What's great about rsync is that it only copies what has changed in the source directory. So if you regularly rsync a large directory to TACC, it may take a long time the 1st time, but the 2nd time (say after downloading more sequencing data to the source), only the new files will be copied.

rsync is a very complicated program, with many options (http://rsync.samba.org/ftp/rsync/rsync.html). However, if you use it like this for directories, it's hard to go wrong:

    rsync -avrP local/path/to/source_directory/ local/path/to/destination_directory/

The -avrP options say "archive mode" (preserve file modification date/time), verbose, recursive and show Progress. Since these are all single-character options, they can be combined after one option prefix dash ( - ).

The trailing slash ( / ) on the source and destination directories are very important! rsync will create the last directory level for you, but earlier levels must already exist.

    rsync -avrP /corral-repl/utexas/BioITeam/web/ucsc_custom_tracks/ data/custom_tracks/

Now repeat the rsync and see the difference:

    rsync -avrP /corral-repl/utexas/BioITeam/web/ucsc_custom_tracks/ $SCRATCH/data/custom_tracks/

#### Copy from a remote computer - scp or rsync

Provided that the remote computer is running Linux and you have ssh access to it, you can use various Linux commands to copy data over a secure connection.

The good news is that once you have learned cp and local rsync, remote secure copy (scp) and remote rsync are very similar!


#### scp

The scp command copies one or more files from a source to a destination, where either source or destination can be a remote path.
Remote paths are similar to local paths, but have user and host information first:

    user_name@full.host.name:/full/path/to/directory/or/file
    -- or –
    user_name@host.name:~/path/relative/to/home/directory

Copy a single file your $SCRATCH/data/test1 directory. We don't really need to access corral remotely, of course, but this shows the remote syntax needed. Be sure to change userid below to your TACC user id!

    scp  userid@stampede.tacc.utexas.edu:/corral-repl/utexas/BioITeam/web/README.txt  $SCRATCH/data/test1
    ls $SCRATCH/data/test1

**NOTE:** The 1st time you access a new host the SSH security prompt will appear, You will always be prompted for your remote host password. The  -r recursive argument works for scp also, just like for cp

#### remote rsync

rsync can be run just like before, but using the remote-host syntax. Here we use two tricks:
-   The tilde ( ~ ) at the start of the path means "relative to my home directory"
-   We traverse through the BioITeam symbolic link created in your home directory earlier.
-   We use the same tilde ( ~ ) in the destination to traverse the scratch symbolic link in your home directory.

Don't forget to change userid below.

    rsync -avrP userid@stampede.tacc.utexas.edu:~/BioITeam/web/ucsc_custom_tracks/ ~/scratch/data/custom_tracks/
    
## EXERCISE

Hit Tab Tab as much as possible to save typing!

    cd
    cp -r /corral-repl/utexas/BioITeam/linuxpractice .
    cd linuxpractice
    cd what
    cat readme

## stamp:~/linuxpractice/what/starts/here/changes/the/world





In this hands-on will learn how to align DNA and RNA-seq data with most widely used software today. Building a whole genome index requires a lot of RAM memory and almost one hour in a typical workstation, for this reason **in this tutorial we will work with chromosome 21** to speed up the exercises. The same steps would be done for a whole genome alignment. Two different datasets, high and low quality have been simulated for DNA, high quality contains 0.1% of mutations and low quality contains 1%. For RNA-seq a 100bp and 150bp datasets have been simulated.


### NGS aligners used:

- [BWA] : BWA is a software package for mapping **DNA** low-divergent sequences against a large reference genome, such as the human genome.
- [HPG Aligner] : HPG Aligner is a new NGS aligner for mapping both **DNA Genomic** and **RNA-seq** data against a large reference genome. It's has been designed for having a high sensitivity and performance.
- [Bowtie2] : *Bowtie 2* is an ultrafast and memory-efficient tool for aligning **DNA** sequencing reads to long reference sequences.
- [TopHat2] : *TopHat* is a fast splice junction mapper for RNA-Seq reads. It aligns **RNA-Seq** reads to mammalian-sized genomes using the ultra high-throughput short read aligner Bowtie, and then analyzes the mapping results to identify splice junctions between exons.
- [STAR] : *STAR* aligns **RNA-seq** reads to a reference genome using uncompressed suffix arrays.

### Other software used in this hands-on:
- [SAMTools] : SAM Tools **provide various utilities** for manipulating alignments in the SAM format, including sorting, merging, indexing and generating alignments in a per-position format.
- [dwgsim] (optional): dwgsim can perform whole **genome simulation**.
- [BEERS] (optional): BEERS is a **simulation engine** for generating **RNA-Seq** data.

### File formats explored:

- [SAM](http://samtools.sourceforge.net/SAMv1.pdf): Sequence alignment format, plain text.
- [BAM](http://www.broadinstitute.org/igv/bam): Binary and compressed version of SAM


### Data used in this practical

Create a ```data``` folder in your **working directory** and download the **reference genome sequence** to be used (human chromosome 21) and *simulated datasets* from **Dropbox** [data](https://www.dropbox.com/sh/4qkqch7gyt888h7/AABD_i9ShwryfAqGeJ0yqqF3a).
For the rest of this tutorial the **working directory** will be **cambridge_mda14** and all the **paths** will be relative to that working directory:
    
    cd cambridge_mda14
    mkdir data

##### Download reference genome from [Ensembl]

Working with NGS data requires a high-end workstations and time for building the reference genome indexes and alignment. During this tutorial we will work only with chromosome 21 to speed up the runtimes. You can download it from **Dropbox** [data](https://www.dropbox.com/sh/4qkqch7gyt888h7/AABD_i9ShwryfAqGeJ0yqqF3a) or from the *Download* link at the top of [Ensembl] website and then to *Download data via FTP*, you get it in only one step by going to:

[http://www.ensembl.org/info/data/ftp/index.html](http://www.ensembl.org/info/data/ftp/index.html)

You should see a species table with a Human (*Homo sapiens*) row and a *DNA (FASTA)* column or click at [ftp://ftp.ensembl.org/pub/release-75/fasta/homo_sapiens/dna/](ftp://ftp.ensembl.org/pub/release-75/fasta/homo_sapiens/dna/), download the chromosome 21 file (*Homo_sapiens.GRCh37.75.dna.chromosome.21.fa.gz*) and move it from your browser download folder to your ```data``` folder:

    mv Homo_sapiens.GRCh37.75.dna.chromosome.21.fa.gz path_to_local_data

**NOTE:** For working with the whole reference genome the file to be downloaded is *Homo_sapiens.GRCh37.75.dna.toplevel.fa.gz*


##### Copy simulated datasets

For this hands-on we are going to use small DNA and RNA-seq datasets simulated from chromosome 21. Data has been already simulated using _dwgsim_ software from SAMtools for DNA and _BEERS_ for RNA-seq. You can copy from the shared resources from **Dropbox** [data](https://www.dropbox.com/sh/4qkqch7gyt888h7/AABD_i9ShwryfAqGeJ0yqqF3a) into your ``data`` directory for this practical session. Preparing the data directory:

    cp path_to_course_materials/alignment/* your_local_data/

Notice that the name of the folders and files describe the dataset, ie. ```dna_chr21_100_hq``` stands for: _DNA_ type of data from _chromosome 21_ with _100_nt read lengths of _high_ quality. Where _hq_ quality means 0.1% mutations and _lq_ quality 1% mutations. Take a few minutes to understand the different files.

**NOTE:** If you want to learn how to simulate DNA and RNA-seq for other conditions go down to the end of this tutorial.

##### Real datasets

For those with access to high-end nodes clusters you can index and simulated whole genome datasets or download real datasets from this sources:
- [1000genomes project](http://www.1000genomes.org/)
- [European Nucleotide Archive (ENA)](https://www.ebi.ac.uk/ena/)
- [Sequence Read Archive (SRA)](http://www.ncbi.nlm.nih.gov/sra)


### Installing SAMtools (Optional, already installed)

Check that it is not installed by executing

    samtools

A list of commands should be printed. If not then proceed with the installation.

Download [SAMtools] from *SF Download Page* link and move to the working directory, then uncompress it.

    mv samtools-0.1.19.tar.bz2 working_directory
    cd working_directory
    tar -jxvf samtools-0.1.19.tar.bz2 
    cd samtools-0.1.19
    make

Check that is correct by executing it with no arguments, the different commands available should be printed. You can also copy it to your ```bin``` folder in your home directory, if bin folder exist, to make it available to the PATH:

    samtools
    cp samtools ~/bin


# Exercise 1: NGS Genomic DNA aligment

In this exercise we'll learn how to download, install, build the reference genome index and align in single-end and paired-end mode with the two most widely DNA aligners: *BWA* and *Bowtie2*. But first, create an ```aligners``` folder to store the software, and an ```alignments``` folder to store the results, create those folders in your *working directory* next to ```data```, you can create both folders by executing:

    mkdir aligners
    mkdir alignments

Now go to ```aligners``` and  ```alignments``` folders and create subfolders for *bwa* and *bowtie* to store the indexes and alignments results:

    cd aligners
    mkdir bwa hpg-aligner bowtie

and

    cd alignments
    mkdir bwa hpg-aligner bowtie
    
**NOTE:** Now your working directory must contain 3 folders: data (with the reference genome of chrom. 21 and simulated datasets), aligners and alignments. Your working directory should be similar to this (notice that aligners have not been downloaded):

```
.
├── aligners
│   ├── bowtie
│   ├── bwa
│   ├── hpg-aligner
├── alignments
│   ├── bowtie
│   ├── bwa
│   ├── hpg-aligner
├── data
│   ├── dna_chr21_100_hq_read1.fastq
│   ├── dna_chr21_100_hq_read2.fastq
│   ├── dna_chr21_100_lq_read1.fastq
│   ├── dna_chr21_100_lq_read2.fastq
│   ├── Homo_sapiens_cDNAs_chr21.fa
│   ├── Homo_sapiens.GRCh37.75.dna.chromosome.21.fa
│   ├── rna_chr21_100_hq_read1.fastq
│   └── rna_chr21_100_hq_read2.fastq
│   ├── rna_chr21_150_lq_read1.fastq
│   └── rna_chr21_150_lq_read2.fastq


```


### BWA
[BWA] is probably the most used aligner for DNA. AS the documentation states it consists of three different algorithms: *BWA*, *BWA-SW* and *BWA-MEM*. The first algorithm, which is the oldest, is designed for Illumina sequence reads up to 100bp, while the rest two for longer sequences. BWA-MEM and BWA-SW share similar features such as long-read support and split alignment, but BWA-MEM, which is the latest, is generally recommended for high-quality queries as it is faster and more accurate. BWA-MEM also has better performance than BWA for 70-100bp Illumina reads.

All these three algorithms come in the same binary so only one download and installation is needed.

##### Download and install (Optional, already installed)

First check that bwa is not currently installed by executing:

    bwa

A list of commands will be printed if already installed. If not you can continue with the installation.

You can click on ```SF download page``` link in the [BWA] page or click directly to:

[http://sourceforge.net/projects/bio-bwa/files](http://sourceforge.net/projects/bio-bwa/files)

Click in the last version of BWA and wait for a few seconds, as the time of this tutorial last version is **bwa-0.7.10.tar.bz2**, the download will start. When downloaded go to your browser download folder and move it to aligners folder, uncompress it and compile it:

    mv bwa-0.7.10.tar.bz2 working_directory/aligners/bwa
    tar -jxvf bwa-0.7.10.tar.bz2
    cd bwa-0.7.10
    make
    cp bwa ~/bin

You can check that everything is allright by executing:

    bwa

Some information about the software and commands should be listed.


##### Build the index

Create a folder inside ```aligners/bwa``` folder called ```index``` to store the BWA index and copy the reference genome into it:
    
    mkdir index
    cp  data/Homo_sapiens.GRCh37.75.dna.chromosome.21.fa  aligners/bwa/index/   (this path can be different!)
    
Now you can create the index by executing:

    bwa index aligners/bwa/index/Homo_sapiens.GRCh37.75.dna.chromosome.21.fa

Some files will be created in the ```index``` folder, those files constitute the index that BWA uses.

 **NOTE:**
 The index must created only once, it will be used for all the different alignments with BWA.


##### Aligning with new BWA-MEM in both single-end (SE) and paired-end (PE) modes

BWA-MEM is the recommended algorithm to use now. You can check the options by executing:

    bwa mem

To align **SE** with BWA-MEM execute:

    bwa mem -t 4 -R "@RG\tID:foo\tSM:bar\tPL:Illumina\tPU:unit1\tLB:lib1" aligners/bwa/index/Homo_sapiens.GRCh37.75.dna.chromosome.21.fa data/dna_chr21_100_hq_read1.fastq > alignments/bwa/dna_chr21_100_hq_se.sam

Now you can use SAMtools to create the BAM file from the *alignment/bwa* folder:
    cd alignments/bwa
    samtools view -S -b dna_chr21_100_hq_se.sam -o dna_chr21_100_hq_se.bam


To align **PE** with BWA-MEM just execute the same command line with the two FASTQ files:

    bwa mem -t 4 -R "@RG\tID:foo\tSM:bar\tPL:Illumina\tPU:unit1\tLB:lib1" aligners/bwa/index/Homo_sapiens.GRCh37.75.dna.chromosome.21.fa data/dna_chr21_100_hq_read1.fastq data/dna_chr21_100_hq_read2.fastq > alignments/bwa/dna_chr21_100_hq_pe.sam
    
Now you can use SAMtools to create the BAM file from the *alignment/bwa* folder:

    cd alignments/bwa
    samtools view -S -b dna_chr21_100_hq_pe.sam -o dna_chr21_100_hq_pe.bam


Now you can do the same for the **low** quality datasets.


##### Aligning with old BWA algorithm, two command lines: ALN and SAMSE/SAMPE in SE and PE modes (Optional exercise)

Now we are going to align SE and PE the **high** quality dataset. Single-end alignment with BWA requires 2 executions. The first uses ```aln``` command and takes the ```fastq``` file and creates a ```sai``` file; the second execution uses ```samse``` and the ```sai``` file and create the ```sam``` file. Results are stored in ```alignments``` folder:

    bwa aln aligners/bwa/index/Homo_sapiens.GRCh37.75.dna.chromosome.21.fa -t 4 data/dna_chr21_100_hq_read1.fastq -f alignments/bwa/dna_chr21_100_hq_se.sai
    bwa samse aligners/bwa/index/Homo_sapiens.GRCh37.75.dna.chromosome.21.fa alignments/bwa/dna_chr21_100_hq_se.sai data/dna_chr21_100_hq_read1.fastq -f alignments/bwa/dna_chr21_100_hq_se.sam


For paired-end alignments with BWA 3 executions are needed: 2 for ```aln``` command and 1 for ```sampe``` command:

    bwa aln aligners/bwa/index/Homo_sapiens.GRCh37.75.dna.chromosome.21.fa -t 4 data/dna_chr21_100_hq_read1.fastq -f alignments/bwa/dna_chr21_100_hq_pe1.sai
    bwa aln aligners/bwa/index/Homo_sapiens.GRCh37.75.dna.chromosome.21.fa -t 4 data/dna_chr21_100_hq_read2.fastq -f alignments/bwa/dna_chr21_100_hq_pe2.sai
    bwa sampe aligners/bwa/index/Homo_sapiens.GRCh37.75.dna.chromosome.21.fa alignments/bwa/dna_chr21_100_hq_pe1.sai alignments/bwa/dna_chr21_100_hq_pe2.sai data/dna_chr21_100_hq_read1.fastq data/dna_chr21_100_hq_read2.fastq -f alignments/bwa/dna_chr21_100_hq_pe.sam

Now you can use SAMtools to create the BAM file from the *alignment/bwa* folder:

    cd alignments/bwa
    samtools view -S -b dna_chr21_100_hq_se.sam -o dna_chr21_100_hq_se.bam
    samtools view -S -b dna_chr21_100_hq_pe.sam -o dna_chr21_100_hq_pe.bam

Now you can do the same for the **low** quality datasets.


### HPG Aligner
[HPG Aligner] is an ultrafast and high sensitivity tool for aligning sequencing reads to long reference sequences. It is particularly good at aligning reads of any size. HPG Aligner 2.0 indexes the genome using a Suffix Arrays. HPG Aligner 2 supports gapped, local, and paired-end alignment modes, together with INDEL realignment and recalibration.

##### Build the index

Create a folder inside HPG Aligner program called ```index``` to store the HPG Aligner index and copy the reference genome into it:
    
    cd hpg-aligner   (_inside aligners folder_)
    mkdir index
    cp data/Homo_sapiens.GRCh37.75.dna.chromosome.21.fa aligners/bwa/index/
    
Now you can create the index by executing:

    hpg-aligner build-sa-index -g aligners/hpg-aligner/index/Homo_sapiens.GRCh37.75.dna.chromosome.21.fa -i aligners/hpg-aligner/index/

Some files will be created in the ```index``` folder, those files constitute the index that HPG Aligner uses.

**NOTE:** The index must created only once, it will be used for all the different alignments with HPG Aligner.


##### Aligning in SE and PE modes


Mapping **SE** with HPG Aligner requires only 1 execution, for aligning the **high** in SE mode execute:

    hpg-aligner dna --cpu-threads 4 -i aligners/hpg-aligner/index/ -f data/dna_chr21_100_hq_read1.fastq -o alignments/hpg-aligner/ --prefix dna_chr21_100_hq_se

And create the BAM file using SAMtools, you could create the BAM file adding _--bam-format_ to the previous command line:

    cd alignments/hpg-aligner
    samtools view -S -b dna_chr21_100_hq_se_out.sam -o dna_chr21_100_hq_se.bam


Mapping in **PE** also requires only one execution:

    hpg-aligner dna --cpu-threads 4 -i aligners/hpg-aligner/index/ -f data/dna_chr21_100_hq_read1.fastq -j data/dna_chr21_100_hq_read2.fastq -o alignments/hpg-aligner --prefix dna_chr21_100_hq_pe

And create the BAM file using SAMtools:

    cd alignments/hpg-aligner
    samtools view -S -b dna_chr21_100_hq_pe_out.sam -o dna_chr21_100_hq_pe.bam
    
Repeat the same steps for the **low** quality dataset.


### Bowtie2

[Bowtie2] as documentation states is an ultrafast and memory-efficient tool for aligning sequencing reads to long reference sequences. It is particularly good at aligning reads of about 50 up to few 100s. Bowtie 2 indexes the genome with an FM Index to keep its memory footprint small: for the human genome, its memory footprint is typically around 3.2 GB. Bowtie 2 supports gapped, local, and paired-end alignment modes.

##### Download and install (Optional, already installed)

First check that bwa is not currently installed by executing:

    bowtie2

A list of commands will be printed if already installed. If not you can continue with the installation.

From [Bowtie2] go to ```Latest Release``` and download the program or go directly to:

[http://sourceforge.net/projects/bowtie-bio/files/bowtie2/2.2.3/](http://sourceforge.net/projects/bowtie-bio/files/bowtie2/2.2.3/)

Click in the Linux version of Bowtie2 and wait for a few seconds, as the time of this tutorial last version is **bowtie2-2.2.3-linux-x86_64.zip**, the download will start. When downloaded go to your browser download folder and move it to aligners folder and uncompress it. No need to compile if you downloaded the Linux version:

    mv bowtie2-2.2.3-linux-x86_64.zip working_directory/aligners/bowtie
    unzip bowtie2-2.2.3-linux-x86_64.zip
    cd bowtie2-2.2.3

You can check that everything is allright by executing:

    bowtie2

Big information about the software and commands should be listed.


##### Build the index

Create a folder inside Bowtie2 program called ```index``` to store the Bowtie2 index and copy the reference genome into it:
    
    cd bowtie   (_inside aligners folder_)
    mkdir index
    cp data/Homo_sapiens.GRCh37.75.dna.chromosome.21.fa aligners/bowtie/index/
    
Now you can create the index by executing:

    bowtie2-build aligners/bowtie/index/Homo_sapiens.GRCh37.75.dna.chromosome.21.fa aligners/bowtie/index/Homo_sapiens.GRCh37.75.dna.chromosome.21.fa

Some files will be created in the ```index``` folder, those files constitute the index that Bowtie2 uses.

**NOTE:** The index must created only once, it will be used for all the different alignments with Bowtie2.


##### Aligning in SE and PE modes

Mapping **SE** with Bowtie2 requires only 1 execution, for aligning the **high** in SE mode execute:

    bowtie2 -q -p 4 -x aligners/bowtie/index/Homo_sapiens.GRCh37.75.dna.chromosome.21.fa -U data/dna_chr21_100_hq_read1.fastq -S alignments/bowtie/dna_chr21_100_hq_se.sam

And create the BAM file using SAMtools;

    cd alignments/bowtie
    samtools view -S -b dna_chr21_100_hq_se.sam -o dna_chr21_100_hq_se.bam


Mapping in **PE** also requires only one execution:

    bowtie2 -q -p 4 -x aligners/bowtie/index/Homo_sapiens.GRCh37.75.dna.chromosome.21.fa -1 data/dna_chr21_100_hq_read1.fastq -2 data/dna_chr21_100_hq_read2.fastq -S alignments/bowtie/dna_chr21_100_hq_pe.sam

And create the BAM file using SAMtools;

    cd alignments/bowtie
    samtools view -S -b dna_chr21_100_hq_pe.sam -o dna_chr21_100_hq_pe.bam
    
Repeat the same steps for the **low** quality dataset.


### More exercises

- Try to simulate datasets with longer reads and more mutations to study which aligner behaves better
- Test the aligner sensitivity to INDELS
- Try BWA-MEM algorithm and compare sensitivity. The same index is valid, only one execution for the SAM file ```./bwa mem index/Homo_sapiens.GRCh37.75.dna.chromosome.21.fa ../../data/dna_chr21_100_low/dna_chr21_100_low.bwa.read1.fastq```


# Exercise 2: NGS RNA-seq aligment

In this exercise we'll learn how to download, install, build the reference genome index and align in single-end and paired-end mode with the two most widely RNA-seq aligner: *TopHat2*. TopHat2 uses Bowtie2 as an aligner.

**NOTE:** Two others commonly used RNA-seq aligners are [STAR] and [MapSplice2], no guided exercises have been documented in this tutorials, but users are encouraged to follow the instructions of their web sites.

Go to ```alignments``` folder and create to folders for *bwa* and *bowtie* to store alignments results:

    cd alignments
    mkdir tophat

**NOTE:** No index is needed for TopHat as it uses Bowtie2 for alignment.
    
### TopHat2

[TopHat2] states to be a *fast* splice junction mapper for RNA-Seq reads, which is not always completrly true. It aligns RNA-Seq reads to mammalian-sized genomes using the ultra high-throughput short read aligner Bowtie, and then analyzes the mapping results to identify splice junctions between exons.

##### Download and install (Optional, already installed)

First check that bwa is not currently installed by executing:

    tophat2

A list of commands will be printed if already installed. If not you can continue with the installation.

From [TopHat2] go to ```Releases``` and download the Linux program by clicking in *Linux x86_64 binary* link.

As the time of this tutorial last version is **tophat-2.0.10.Linux_x86_64.tar.gz**, the download will start. When downloaded go to your browser download folder and move it to aligners folder and uncompress it. No need to compile if you downloaded the Linux version:

    mv tophat-2.0.10.Linux_x86_64.tar.gz working_directory/aligners/tophat
    tar -zxvf tophat-2.0.10.Linux_x86_64.tar.gz
    cd tophat-2.0.10.Linux_x86_64

You can check that everything is allright by executing:

    tophat2

Big information about the software and commands should be listed.

**NOTE:** TopHat uses Bowtie as the read aligner. You can use either Bowtie 2 (the default) or Bowtie (--bowtie1) and you will need the following Bowtie 2 (or Bowtie) programs in your PATH. Index must be created with Bowtie not TopHat. So, copy Bowtie2 into ~/bin:
    
    cd bowtie2   (bowtie 2.2 does not work)
    cp bowtie* ~/bin


##### Aligning in SE and PE modes

To align in SE mode:

    tophat2 -o alignments/tophat/rna_chr21_100_hq_se aligners/bowtie/index/Homo_sapiens.GRCh37.75.dna.chromosome.21.fa data/rna_chr21_100_hq_read1.fastq

And for PE:

    tophat2 -o alignments/tophat/rna_chr21_100_hq_pe/ aligners/bowtie/index/Homo_sapiens.GRCh37.75.dna.chromosome.21.fa data/rna_chr21_100_hq_read1.fastq data/rna_chr21_100_hq_read2.fastq

Now align the rna dataset of 150bp with low quality and compare stats.


### STAR and MapSplice2
[STAR] and [MapSplice2] are two others interesting RNA-seq aligners. [STAR] offer a great performance while still have good sensitivity. [MapSplice2] shows usually a better sensitivity but is several times slower.

##### STAR installation (Optional, already installed)
STAR comes compiled for Linux, you only have to download the *tarball*:

    tar -zxvf STAR_2.3.0e.Linux_x86_64_static.tgz

Read the documentation and try to align the simulated dataset.


##### MapSplice2 installation (Optional, already installed)
MapSplice must be unizpped and compiled:

    unzip MapSplice-v2.1.6.zip
    cd MapSplice-v2.1.6
    make

Read the documentation and try to align the simulated dataset.


# Exercise 3: Simulating NGS datasets (Optional)

### DNA
Download [dwgsim] from http://sourceforge.net/projects/dnaa/files/ to the *working_directory* and uncompress it and compile it:

    tar -zxvf dwgsim-0.1.10.tar.gz
    cd dwgsim-0.1.10
    make

Check options by executing:

    ./dwgsim

Then you can simulate 2 million reads of 150bp with a 2% if mutation executing:

    ./dwgsim-0.1.11/dwgsim -1 150 -2 150 -y 0 -N 2000000 -r 0.02 ../data/Homo_sapiens.GRCh37.75.dna.chromosome.21.fa ../data/dna_chr21_100_low/dna_chr21_100_verylow


### RNA-seq
[BEERS] is a perl-based program, no compilation is needed, just download it from here http://www.cbil.upenn.edu/BEERS and uncompress it:

    tar xvf beers.tar

