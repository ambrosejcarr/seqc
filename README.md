## SEquence Quality Control (SEQC -- /sek-si:/)

### Overview:

SEQC is a package that is designed to process sequencing data on the cloud. It also
contains tools to analyze processed data, which has a much smaller memory footprint. 
Thus, it should be installed both locally and on your remote server. However, it should
be noted that some portions of the data pre-processing software require 30GB of RAM to 
run, and are unlikely to work on most laptops and workstations.
 
To faciliate easy installation and use, we have made available Amazon Machine Images
(AMIs) that come with all of SEQC's dependencies pre-installed. In addition, we have
uploaded the indices (-i/--index parameter, see "Running SEQC) and barcode data
(-b/--barcodes) to public amazon s3 repositories. These links can be provided to SEQC and
it will automatically fetch them prior to initiating an analysis run. 

Amazon Web Services (AWS) is only accessible through a programmatic interface. To
simplify the starting and stopping of amazon compute servers, we have written several
plug-ins for a popular AWS interface, StarCluster. By installing starcluster, users can
simply call starcluster start <cluster name>; starcluster sm <cluster name> for instant
access to a compute server that is immediately ready to run SEQC on data of your choosing.

### Installation \& Dependencies:

#### Dependencies For Remotely Running on AWS:
1. Amazon Web Services (optional): If you wish to run SEQC on amazon (AWS), 
you must set up an AWS account (see below). SEQC will function on any machine running
a nix-based operating system, but we only provide AMIs for AWS. 
3. Starcluster (optional; install locally): If you run SEQC on aws, then you create aws
compute instances with the python 2.7 
<a href=https://github.com/jtriley/StarCluster>starcluster</a>
package, which streamlines creation and shutdown of servers that are ready to run SEQC.
Requires <a href=https://www.python.org/downloads/release/python-2710/>Python 2.7</a> 
with <a href=http://pip.readthedocs.org/en/stable/installing/>pip</a>.
            
        $> git clone git://github.com/jtriley/StarCluster.git
        $> cd StarCluster
        $> sudo python distribute_setup.py
        $> sudo python setup.py install
    
3. Install Starcluster plugins and config template from SEQC:

        $> git clone https://github.com/ambrosejcarr/seqc.git
        $> cp seqc/src/plugins/*.py ~/.starcluster/plugins/
        $> cp seqc/src/plugins/starcluster.config ~/.starcluster/config

#### Dependencies for Local Installation or Other Cloud Computing Platforms:
1. <a href=https://www.python.org/downloads/>Python 3</a>
2. <a href=https://www.hdfgroup.org/HDF5>libhdf5</a>, a highly efficient database used to
store SEQC output.
3. SEQC depends on several python3 packages that are automatically installed and updated.
to view these packages, please view the `setup.py` file packaged with SEQC.

### Setting up HDF5 on your local computer:
1. After downloading libhdf5 from source, it can be installed by typing:
```
    $> ./configure --prefix=/usr/local/
    $> make
    $> make install
```
2. Then, you can install pytables by typing:
    ```pip3 install tables```
3. If you installed libhdf5 without giving arguments in the "configure" step, you can also:
    ```export HDF_DIR=/your/installation/directory/for/hdf5```
But before you do that, you need to make sure that you have the prereqs installed previously:
    * numpy
    * numexpr
    * cython

### Setting up AWS, SEQC, and starcluster

Once all dependencies have been installed, SEQC can be installed on any machine by typing:

    $> git clone https://github.com/ambrosejcarr/seqc.git
    $> pip3 install -e seqc/

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

#### Personalize the Dummy StarCluster Config File Provided by SEQC.
1. Open the `~/.starcluster/config` file
2. Under `[aws info]` enter the following information: 
    1. `AWS_ACCESS_KEY_ID = #access_id` (This is your AWS Access ID from Step (1))
    2. `AWS_SECRET_ACCESS_KEY = #secret_key` (This is your AWS Secret Key from Step (1))
    3. `AWS_USER_ID= #user_id` (This is a numerical ID from AWS, found under IAM users)
    4. Click on your username on the top right corner of the AWS dashboard and click
    “My Account” -- your Account Id should pop up at the top of the page (a 12-digit 
    number)
3. Under Defining EC2 Keypairs:
    1. rename `[key <your_key_name>]` to the name of the key you generate above.
    2. change key location to the location of your `<keyname.rsa>` file:
    `KEY_LOCATION=~/.ssh/<keyname>.rsa`
5. Under templates, find `[cluster c3.large]`
    1. change key to `<your_key_name>`
    
#### Install and Configure AWS CLI (AWS Command Line Interface).
1. You can install by typing `pip install awscli`
2. Then, configure it by typing `aws configure`:
    * AWS Access Key ID [*******]: `access_id`
    * AWS Secret Access Key [*******]: `secret_key`
    * Default region name [us-west-2]: `us-east-1` (Adjust accordingly)
    * Default output format [None]: `text`

#### Start a cluster:
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


### Running SEQC:

After SEQC is installed, help can be listed:

    $> SEQC -h
    usage: SEQC [-h]
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
