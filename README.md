# SEquence Quality Control

## Dependencies:

Most dependencies will be automatically installed by SEQC. However, there are a few
external dependencies that must be installed in order for certain functions to be enabled.

1. Amazon Web Services (optional): If you wish to run SEQC on amazon (AWS), 
you must set up an AWS account (see below). 
2. Starcluster (optional): If you run SEQC on aws, that you create aws compute instances 
with the python 2.7 <a href=https://github.com/jtriley/StarCluster>starcluster</a>
package, which streamlines creation and shutdown of servers that are ready to run SEQC.
2. libhdf5: We have provided amazon machine images or "AMIs" with all of starcluster's 
external dependencies pre-installed. However, if you wish to use your own machine image, 
you will need to install <a href=https://www.hdfgroup.org/HDF5>libhdf5</a>. Similarly, if
you wish to parse any intermediate data offline on your local machine, you will need to 
install libhdf5. 
3. SEQC depends on several python3 packages that are automatically installed and updated.
to view these packages, please view the `setup.py` file packaged with SEQC.

## Setting up SEQC and starcluster

### Creating a new AWS account:
1. Navigate <a href=http://aws.amazon.com>here</a> and click “Create an AWS Account.”
2. Enter your desired login e-mail and click the “I am a new user” radio button.
3. Fill out the Login Credentials with your name and password.
4. Fill out your contact information, read the AWS Customer Agreement, and accept if you
wish to create an account.
5. Save your AWS Access ID and Secret Key -- this information is very important!

### Create an RSA key to launch a new cluster:
1. Sign into your AWS account and go to the EC2 Dashboard.
2. Click “Key Pairs” in the NETWORK & SECURITY tab.
3. Click “Create Key Pair” and give it a new name.
4. This will install a new key called <keyname>.pem on your local machine. 
5. Rename the key to an .rsa extension and move it to a secure directory.
6. example: `<keyname>.rsa` and move it to a directory (e.g. `~/.ssh/`)
7. Change the permission settings with `$> chmod 600 /path/to/key/keyname.rsa`

### Next, install StarCluster by cloning the Github repo.
1. `$> git clone https://github.com/ambrosejcarr/StarCluster.git`
2. `$> sudo python setup.py install`

### Set up the default StarCluster config file.
1. Type `starcluster help` and enter `2`
2. Under `[aws info]` enter the following information: 
    1. `AWS_ACCESS_KEY_ID = #access_id` (This is your AWS Access ID from Step (1))
    2. `AWS_SECRET_ACCESS_KEY = #secret_key` (This is your AWS Secret Key from Step (1))
    3. `AWS_USER_ID= #user_id`
    4. Click on your username on the top right corner of the AWS dashboard and click
    “My Account” -- your Account Id should pop up at the top of the page (a 12-digit 
    number)
    5. `AWS_REGION_NAME = us-east-1a` (or any other of 1b, 1c, 1d, 1e - they have 
    different spot prices) 
3. Under Defining EC2 Keypairs:
    1. Create a section called: `[key <keyname>]`
    2. below, write: `KEY_LOCATION=~/.ssh/<keyname>.rsa`
    3. Under `[cluster smallcluster]`:
    4. `KEYNAME = <keyname>`
    5. `NODE_IMAGE_ID = %(x86_64_ami)s`
    6. `SUBNET_ID=subnet-99999999`
        1. Go to the “VPC” under the Networking Tab in the AWS Dashboard
        2. Click on “Subnets”
        3. Find the appropriate subnet according to your availability zone
        (e.g. “subnet-99999999”) and write it there
    7. `AVAILABILITY_ZONE = us-east-1a` (or other, if you changed this in 2.5)

### Start a cluster:
1. `$> starcluster start -c <template_name> <cluster_name>`
2. Wait until the cluster is finished setting up. Then, the cluster can be accessed
using:
3. `$> starcluster sshmaster -X <cluster_name>`  (-X gives you x-window plotting capability)
4. To exit the cluster, simply type “exit”.
5. Other things like `starcluster stop <cluster_name>`, `terminate`, `start -c`, etc.
6. You can also copy files to/from the cluster using the put and get commands. 
To copy a file or entire directory from your local computer to the cluster:
`$> starcluster put mycluster /path/to/local/file/or/dir /remote/path/`
7. To copy a file or an entire directory from the cluster to your local computer:
`$> starcluster get mycluster /path/to/remote/file/or/dir /local/path/`


## Running SEQC:

After SEQC is installed, help can be listed:

    $> SEQC -h
    usage: SEQC [-h] {in-drop,drop-seq,mars-seq,cel-seq,avo-seq,strt-seq} ...
    
    positional arguments:
      {in-drop,drop-seq,mars-seq,cel-seq,avo-seq,strt-seq}
                            library construction method types
        in-drop             in-drop help
        drop-seq            drop-seq help
        mars-seq            mars-seq help
        cel-seq             cel-seq help
        avo-seq             avo-seq help
        strt-seq            strt-seq help
    
    optional arguments:
      -h, --help            show this help message and exit


Help on parsing individual data types can be obtained by typing:

    $> SEQC in-drop -h
    usage: SEQC in-drop [-h] [-i I] [-n N] [-o O] [-b B] [-f [F [F ...]]]
                        [-r [R [R ...]]] [-s [S]] [-m M] [-l L]
                        [--star-args SA [SA ...]] [--list-default-star-args]
    
    optional arguments:
      -h, --help            show this help message and exit
    
    Required Arguments:
      -i I, --index I       star alignment index folder. This folder will be
                            created if it does not exist
      -n N, --n-threads N   number of threads to run
      -o O, --output-file O
                            stem of filename in which to store output
      -b B, --barcodes B    location of serialized barcode object.
    
    Input Files:
      pass one input file type: sam (-s), raw fastq (-f, [-r]), or processed
      fastq (-m)
    
      -f [F [F ...]], --forward [F [F ...]]
                            forward fastq file(s)
      -r [R [R ...]], --reverse [R [R ...]]
                            reverse fastq file(s)
      -s [S], --sam [S]     sam file(s) containing aligned, pre-processed reads
      -m M, --merged-fastq M
                            fastq file containing merged, pre-processed records
    
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
observing a co-alignment to other genes.

Human and mouse indices can be found on our aws s3 bucket at
`s3://dplab-data/genomes/mm38/` and `s3://dplab-data/genomes/hg38`. These indices
are built from recent releases of ENSEMBL genomes.

If new indices must be generated, these can be produced by calling:

    $> INDEX=<index_folder>
    $> # valid organisms are hg38 and mm38. Can also produce joint indicies with 
    $> # "[hg38, mm38]"
    $> ORG=<organism>
    $> # STAR will use this number of threads to produce the index
    $> THREADS=<number of threads to use>  
    $> mkdir $INDEX  # create index directory
    $> python3 -c "import seqc.align; seqc.align.STAR.build_index($INDEX, $ORG, $THREADS)"

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
