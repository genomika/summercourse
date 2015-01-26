% [WorkShop de Verão Bioinformática](http://github.com/genomika/summercourse/)
% __Introdução ao Linux__
% _(updated 25-01-2014)_

<!-- COMMON LINKS HERE -->

# Acessando servidor remoto

## A Janela de Terminal
- Macs e Linux possuem programas de terminal nativos - encontre-os em sua máquina
- Usuários Windows vão precisar de ajuda, opções a seguir:
    - [Putty](http://www.chiark.greenend.org.uk/~sgtatham/putty/download.html) "Putty"
    - [Git-bash](http://msysgit.github.io/) "GitBash"
    - [Cygwin](http://www.cygwin.com/) "Cygwin"

## SSH

ssh é um programa executável que roda em seu computador local e permite que você se conecte de forma segura a um computador remoto. Em Macs, Linux e Windows (Git-Bash out Cygwin), você pode rodar de uma tela de terminal. Responda yes a questão de segurança do prompt do SSH.

##### Acessando servidor remoto via SSH.
    ssh seuusuario@172.16.225.177

- Se você está usando Putty como seu terminal no Windows:
    -  Clique Duplo no icone Putty.exe 
    -  Na janela de configuração do PuTTy 
    -  Garanta que o tipo de conexão é SSH
    -  Digite o ip do servidor 172.16.225.177 (host name)
    -  Clique no botão Open 
    -  Responda  Yes na pergunta de segurança do SSh
    -  No terminal putty
    -  Digite seu nome de usuario quando solicitado e enter


## Shell Bash

Você está agora na linha de comando! É como se você estivesse executando comandos direto no computador remoto, mas na verdade há 2 programas se comunicando: seu terminal local e o shell remoto. Há vários programas shell disponíveis em linux, mas o padrão é o bash (Bourne-again shell). O Terminal é bastante "básico" – apenas envia o que você digita através de uma conexão por SSL (secure sockets layer) ao servidor, e finalmente exibindo o texto de volta ao shell.  O trabalho real é feito no computador remoto, por programas chamados no shell bash.

![Terminal Example](https://wikis.utexas.edu/download/attachments/66696867/terminal.png?version=1&modificationDate=1400007439000&api=v2)

##### Preparando seu ambiente de trabalho.

Primeiro crie alguns diretórios e links símbolicos que nós iremos usar (falaremos sobre isto depois).

**NOTE:** Você pode copiar e colar estas linhas abaixo de código na janela do seu terminal. Apenas lembre de pressionar 'Enter' após a última linha.

    cd
    mkdir work
    ln -s -f /usr/local/share/tools/ tools
    ln -s -f /mnt/data/ data

- The mkdir command creates a new directory.
- The ln -s command creates a symbolic link, a shortcut the the linked file or directory.
- here the link targets are your work file system areas
- having this link shortcut will help when you want to copy files to your work, and when you navigate the file system using a remote SFTP client
- always change directory (cd) to the directory where we want the links created before executing ln -s 
- here we want the links under your home directory (cd with no arguments)

##### Want to know where a link points to? Use ls with the -l (long listing) option.
    ls -l

Set up a $HOME/local/bin directory and link some scripts there that we will use a lot in the class.

    mkdir -p $HOME/local/bin
    cd $HOME/local/bin
    ln -s -f /usr/local/share/tools/bin/cutadapt
    ln -s -f /usr/local/share/tools/bin/samstat

-  The mkdir command creates a new directory. The -p option says to create intermediate directories if needed (like local here).
-  here we're creating a $HOME/local/bin directory where we'll put some programs used in the course
-  $HOME is an environment variable set by Unix that refers to your home directory.
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
    cp /usr/local/common/ngs_profile .profile
    chmod 600 .profile


-  The chmod 600 .profile_user command marks the file as readable and writable only by you. The .profile script file will not be executed unless it has these exact permissions settings.
-  The well-known filename is .profile (or .profile on some systems), which is specific to the bash shell.

Since .profile is executed when you login, to ensure it is set up properly you should first log off  like this:

    exit
    
Then log back in to 172.16.225.177. This time your .profile_user will be executed and you should see a new shell prompt:
    
    ~$

The great thing about this prompt is that it always tells you where you are, which avoids having to issue the pwd (present working directory) command all the time. Execute these commands to see how the prompt reflects your current directory. (Don't just copy-and-paste here because we've included the prompt.)

    $ mkdir -p tmp/a/b/c
    $ cd tmp/a/b/c
    $ ~/tmp/a/b/c$

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

    cat .profile
    

**NOTE:**   The cat command just echos the entire file's content, line by line, without pausing, so should not be used to display large files. Instead, use a "pager" (like more or less) or look at parts of the file with **head** or **tail**.

You'll see the following (you may need to scroll up a bit to see the beginning):

    #!/bin/basH
    # Change the command line prompt to contain the current directory path
    PS1='workshop:\w$ '
    # Try to ensure all created files can be read/writtin by group members
    umask 002
    # Add current directory and $HOME/local/bin to PATH
    export PATH=.:$HOME/local/bin:$PATH
    # Use yellow for directories, not that horrible blue
    dircolors .dircolors > /dev/null


So what does the common profile file do? Several things. Let's look at a few of them.

#### she-bang

The first line is the "she-bang". It tells the shell what program should execute this file – in this case, bash itself – even though the expression is inside a shell comment (denoted by the # character).

    #!/bin/bash

#### environment variables

    # Environment variables for useful locations
    export PATH=.:$HOME/local/bin:$PATH

Environment variables are like variables in a programming language like python or perl (in fact bash is a complete programming language). They have a name (like PATH above) and a value (the value for PATH is the pathname **$HOME/local/bin:$PATH**).

#### shell completion

You can use these environment variables to shorten typing, for example, to look at the contents of the shared  directory as shown below, using the magic Tab key to perform shell completion.

    # now hit Tab twice to see the contents of the directory
    ls /genomika/bioinformatics/
    # now type "co" and hit Tab again
    ls /genomika/bioinformatics/co
    # your command line should now look like this
    ls /genomika/bioinformatics/core_ngs_tools/
    # now type "m" and one Tab
    ls /genomika/bioinformatics/core_ngs_tools/m
    # now just type one Tab
    ls /genomika/bioinformatics/core_ngs_tools/misc/
    # the shell expands as far as it can unambiguously, so your command line should look like thi
    ls  /genomika/bioinformatics/core_ngs_tools/misc/smallExample
    # type a period (".") then hit Tab twice again -- you're narrowing down the choices
    ls  /genomika/bioinformatics/core_ngs_tools/misc/smallExample
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

As you can see, there are a lot of locations on the **$PATH**. That's because when you load modules at the server (such as the module load lines in the common profile), that mechanism makes the programs available to you by putting their installation directories on your $PATH. We'll learn more about modules shortly.

Here's how the shared profile adds your **$HOME/local/bin** directory to the location list – recall that's where we linked some programs we'll use – along with a special dot character ( . ) that means "here", or "whatever the current directory is".

    # Add current directory and $HOME/local/bin to PATH 
    export PATH=.:$HOME/local/bin:$PATH


### Download from a link – wget

Well, you don't have a desktop at the server to "Save as" to, so what to do with a link? The wget program knows how to access web URLs such as http, https and ftp.

#### wget

Get ready to run wget from the directory where you want to put the data. Don't press Enter after the wget command – just put a space after it.

    cd $WORK
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

    mkdir -p $work/data/test1
    cp  //172.16.225.20/ngscourse/data  $work/data/test1/
    ls $work/data/test1

Copy a directory to your work area. The -r argument says "recursive".

    cds
    cd data
    cp -r  //172.16.225.20/ngscourse/data/general/ general/


#### local rsync

The rsync command is typically used to copy whole directories. What's great about rsync is that it only copies what has changed in the source directory. So if you regularly rsync a large directory to the server, it may take a long time the 1st time, but the 2nd time (say after downloading more sequencing data to the source), only the new files will be copied.

rsync is a very complicated program, with many options (http://rsync.samba.org/ftp/rsync/rsync.html). However, if you use it like this for directories, it's hard to go wrong:

    rsync -avrP local/path/to/source_directory/ local/path/to/destination_directory/

The -avrP options say "archive mode" (preserve file modification date/time), verbose, recursive and show Progress. Since these are all single-character options, they can be combined after one option prefix dash ( - ).

The trailing slash ( / ) on the source and destination directories are very important! rsync will create the last directory level for you, but earlier levels must already exist.

    rsync -avrP /corral-repl/utexas/BioITeam/web/ucsc_custom_tracks/ data/custom_tracks/

Now repeat the rsync and see the difference:

    rsync -avrP /corral-repl/utexas/BioITeam/web/ucsc_custom_tracks/ $work/data/custom_tracks/

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
