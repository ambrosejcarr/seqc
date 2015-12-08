## SEquence Quality Control (SEQC -- /sek-si:/)

### Overview:

SEQC is a package that is designed to process sequencing data on the cloud. It also
contains tools to locally analyze processed data, which has a much smaller memory
footprint. We have created a pre-formatted machine image on Amazon's Web Services (AWS)
with the current SEQC version installed. The user can install and call SEQC on their
local machine, which will trigger an AWS server to start up, process the specified
sequencing experiment, and email you the results when it completes. Please note that the
use of this capability requires the user have an active AWS account. If Users wish to
run SEQC on their own server, they can install SEQC locally and run it without the
--remote flag.

### Requirements:
1. 32GB RAM for human genome processing. Other larger genomes may require more RAM.
2. 4+ compute cores (30+ optimal)

### Installation \& Dependencies:

#### Installation for running remotely:
1. Create an Amazon Web Services account (see below).
2. Create an RSA key and register it with AWS (see below). 
3. Download and install <a href=https://www.python.org/downloads/>Python 3</a>.
4. Install <a href=https://www.hdfgroup.org/HDF5>libhdf5</a> (optional), a highly
efficient database used to store SEQC output. This is only required for remote-running if
you wish to parse the filtered data.
5. Install SEQC:

        $> git clone https://github.com/ambrosejcarr/seqc.git
        $> pip3 install -e seqc/

### Installing HDF5:
#### Installing from Source:
1. After downloading libhdf5 source, it can be installed as follows:

        $> ./configure --prefix=/usr/local/
        $> make
        $> make install

2. Install pytables by typing: `pip3 install tables`

#### Installing without previous configuration
1. If you installed libhdf5 without giving arguments in the "configure" step, make sure
that you have the necessary prereqs already installed:
    * numpy
    * numexpr
    * cython
2. Then set the $HDF_DIR environment variable by typing: 

        $> export HDF_DIR=/your/installation/directory/for/hdf5

3. You should now be able to install pytables: `pip3 install tables`

### Miscellaneous set-up and installation:

#### Setting up an AWS Account: 
1. Navigate <a href=http://aws.amazon.com>here</a> and click “Create an AWS Account.”
2. Enter your desired login e-mail and click the “I am a new user” radio button.
3. Fill out the Login Credentials with your name and password.
4. Fill out your contact information, read the AWS Customer Agreement, and accept if you
wish to create an account.
5. Save your AWS Access ID and Secret Key -- this information is very important!

#### Create an RSA key to allow you to launch a cluster
1. Sign into your AWS account and go to the EC2 Dashboard.
2. Click “Key Pairs” in the NETWORK & SECURITY tab.
3. Click “Create Key Pair” and give it a new name.
4. This will install a new key called <keyname>.pem on your local machine. 
5. Rename the key to an .rsa extension and move it to a secure directory.
6. example: `<keyname>.rsa` and move it to a directory (e.g. `~/.ssh/`)
7. Change the permission settings with `$> chmod 600 /path/to/key/keyname.rsa`

### Running SEQC:

After SEQC is installed, help can be listed:

    $usage: SEQC [-h] [-v]
                {in-drop,drop-seq,mars-seq,cel-seq,avo-seq,strt-seq,index} ...
    
    positional arguments:
      {in-drop,drop-seq,mars-seq,cel-seq,avo-seq,strt-seq,index}
                            library construction method types
        in-drop             in-drop help
        drop-seq            drop-seq help
        mars-seq            mars-seq help
        cel-seq             cel-seq help
        avo-seq             avo-seq help
        strt-seq            strt-seq help
        index               SEQC index functions
    
    optional arguments:
      -h, --help            show this help message and exit
      -v, --version         print version and exit

Help on running each data types can be obtained by typing:

    $> SEQC in-drop -h
    usage: SEQC in-drop [-h] [-i I] [-n N] [-o O] [-c C] [-k K] [--remote]
                        [--email-status E] [--no-terminate] [-b B]
                        [-f [F [F ...]]] [-r [R [R ...]]] [-s [S]] [-m M]
                        [--basespace BS BS] [-l L] [--star-args SA [SA ...]]
                        [--list-default-star-args]
    
    optional arguments:
      -h, --help            show this help message and exit
    
    Required Arguments:
      -i I, --index I       local or s3 location of star alignment index folder.
                            This folder will be created if it does not exist
      -n N, --n-threads N   number of threads to run
      -o O, --output-prefix O
                            stem of filename in which to store output
      -c C, --cluster-name C
                            optional name for aws cluster
      -k K, --aws-upload-key K
                            upload processed results to this AWS folder. Required
                            if --remote is passed
      --remote              run the requested SEQC command remotely
      --email-status E      email results to this address
      --no-terminate        Decide if the cluster will terminate upon completion.
      -b B, --barcodes B    local or s3 location of serialized barcode object.
    
    Input Files:
      pass one input file type: sam (-s), raw fastq (-f, [-r]), or processed
      fastq (-m)
    
      -f [F [F ...]], --forward [F [F ...]]
                            s3 link or filesystem location of forward fastq
                            file(s)
      -r [R [R ...]], --reverse [R [R ...]]
                            s3 link or filesystem location of reverse fastq
                            file(s)
      -s [S], --samfile [S]
                            s3 link or filesystem location of sam file(s)
                            containing aligned, pre-processed reads
      -m M, --merged-fastq M
                            s3 link or filesystem location of fastq file
                            containing merged, pre-processed records
      --basespace BS BS     takes 2 arguments: a basespace sample id, and a
                            basespace oAuth token. If provided, fastq input data
                            will be downloaded from BaseSpace for processing.
    
    Optional arguments for disambiguation:
      -l L, --frag-len L    the number of bases from the 3 prime end to consider
                            when determining trancript overlaps
    
    Optional arguments for STAR aligner:
      --star-args SA [SA ...]
                            additional arguments for STAR. Pass as arg=value
                            without leading "--". e.g. runMode=alignReads
      --list-default-star-args
                            list SEQDB default args for the STAR aligner


All SEQC runs require that you pass a SEQC index (`-i/--index`). These are STAR indices,
augmented by SEQC-specific files:

1. `annotations.gtf`: a modified GTF file, containing truncated sequence sizes that
reflect that we expect all data to fall within ~ 1kb of the transcriptional termination
sites. In addition, transcripts are tagged with "SCIDs", identifiers that merge
transcripts and genes which cannot be distinguished in the ~ 1kb of sequence that we
expect to observe.
2. `p_coalignments_array.p`: a binary file containing, for each SCID, the probability of
observing a co-alignment to other genes or transcripts.

Human and mouse indices can be found on our aws s3 bucket at
`s3://dplab-data/genomes/mm38/` and `s3://dplab-data/genomes/hg38`. These indices
are built from recent releases of ENSEMBL genomes. These links can be passed directly to
SEQC, which will download them before beginning the analysis

If new indices must be generated, these can be produced by the SEQC index method:

    $> SEQC index -h
    usage: SEQC index [-h] [-b] [-t] -o O [O ...] -i I [-n N] [--phix]
    
    optional arguments:
      -h, --help            show this help message and exit
      -b, --build           build a SEQC index
      -t, --test            test a SEQC index
      -o O [O ...], --organism O [O ...]
                            build index for these organism(s)
      -i I, --index I       name of folder where index should be built or
                            containing the index to be verified
      -n N, --n-threads N   number of threads to use when building index
      --phix                add phiX to the genome index and GTF file.
     
     $> # for example, to build a mouse index with phiX features added to mm38, in a
     $> $ folder called 'mouse', using 7 threads
     $> SEQC index -b -o mm38 -i mouse -n 7 --phix
    
Some data types require serialized barcode objects (`-b/--barcodes`). These objects contain
all of the barcodes for an experiment, as they would be expected to be observed.
For example, if you expect to observe the reverse complement of the barcodes you used to
construct the library, then this object should be built from reverse complements.   
 
These barcode files can be found at `s3://dplab-data/barcodes/`. If you need to generate
a new barcode object, this can be accomplished with the built-in `PROCESS_BARCODES`
utility:

    $> PROCESS_BARCODES -h
    usage: PROCESS_BARCODES [-h] [-o O] [-b B [B ...]] [-p P]
                            [--reverse-complement]
    
    optional arguments:
      -h, --help            show this help message and exit
      -o O, --output_stem O
                            name and path for output file
      -b B [B ...], --barcode_files B [B ...]
                            barcode files
      -p P, --processor P   type of experiment barcodes will be used in
      --reverse-complement  indicates that barcodes in fastq files are reverse
                            complements of the barcodes found in barcode files

Example usage:

`$> PROCESS_BARCODES -o ./in_drop_barcodes -b <barcode_file> -p in-drop --reverse-complement`
would save a new, reverse-complemented barcode object built from `<barcode_file>` at
`./in_drop_barcodes.p`

### Inputting Data

SEQC can run on many different types of data, from many different locations. It supports
running from:

1. Raw, unprocessed, paired-end fastq files. These can be passed as individual files using
the `-f` and `-r` parameters, AWS S3 links to individual files, or AWS S3 links to 
folders containing multiple fastq files. In addition, passing a BaseSpace sample ID to 
the `--basespace` command will download your forward and reverse fastq from your BaseSpace 
account. 

2. Merged fastq (`-m`). These files have been demultiplexed and the barcode IDs , UMIs, 
and quality information of the read holding barcode data have been prepended to the read
holding genomic data

3. Sam files (`-s`). These files must be in the proper demultiplexed format (see Merged  
fastq above)

