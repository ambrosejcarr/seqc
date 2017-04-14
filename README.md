## SEquence Quality Control (SEQC -- /sek-si:/)

## Overview:

SEQC is a python package that processes single-cell sequencing data in the cloud and analyzes it interactively on your local machine.
 
To faciliate easy installation and use, we have made available Amazon Machine Images (AMIs) that come with all of SEQC's dependencies pre-installed. In addition, we have uploaded common genome indices (-i/--index parameter) and barcode data (--barcode-files) to public amazon s3 repositories. These links can be provided to SEQC and it will automatically fetch them prior to initiating an analysis run. Finally, it can fetch input data directly from BaseSpace or amazon s3 for analysis.

For users with access to in-house compute clusters, SEQC can be installed on your systems and run using the --local parameter.

### Dependencies:


#### Python3
Python must be installed on your local machine to run SEQC. We recommend installing python3 through your unix operating system's package manager. For Mac OSX users we recommend <a href=http://brew.sh/>homebrew</a>. Typical installation commands would be: 

    brew install python3  # mac
    apt-get install python3  # debian
    yum install python3 # rpm-based


#### STAR
To run SEQC locally, you must install the <a href=https://github.com/alexdobin/STAR>STAR Aligner</a>, <a href=http://www.htslib.org/>Samtools</a>, and <a href=https://support.hdfgroup.org/HDF5/>hdf5</a>. If you only intend to run remotely, these dependency is optional. We recommend installing samtools and hdf5 through your package manager, if possible. 


#### Hardware Requirements:
For processing a single lane (~200M reads) against human- and mouse-scale genomes, SEQC requires 30GB RAM, approximately 200GB free hard drive space, and scales linearly with additional compute cores. If running on AWS (see below), jobs are automatically scaled up or down according to the size of the input. There are no hardware requirements for the computer used to launch remote instances.  


#### Amazon Web Services:
SEQC can be run on any unix-based operating system, however it also features the ability to automatically spawn Amazon Web Services instances to process your data. If you wish to take advantage of AWS, you will need to follow their instructions to:

1. <a href=http://aws.amazon.com>Set up an AWS account</a>
2. <a href=https://aws.amazon.com/cli/>Install and configure AWS CLI</a> 
3. <a href=http://docs.aws.amazon.com/AWSEC2/latest/UserGuide/ec2-key-pairs.html>Create and upload an rsa-key for AWS</a> 


### SEQC Installation:

Once all dependencies have been installed, SEQC can be installed on any machine by typing:

    $> git clone https://github.com/ambrosejcarr/seqc.git
    $> cd seqc && python3 setup.py install 

Please note that to avoid passing the -k/--rsa-key command when you execute SEQC runs, you can also set the environment variable `AWS_RSA_KEY` to the path to your newly created key.

### Running SEQC:

After SEQC is installed, help can be listed:

    SEQC [-h] [-v] {run,progress,terminate,instances,start,index} ...
    
    Processing Tools for scRNA-seq Experiments
    
    positional arguments:
      {run,progress,terminate,instances,start,index}
        run                 initiate SEQC runs
        progress            check SEQC run progress
        terminate           terminate SEQC runs
        instances           list all running instances
        start               initialize a seqc-ready instance
        index               create a SEQC index
    
    optional arguments:
      -h, --help            show this help message and exit
      -v, --version         show program's version number and exit

In addition to processing sequencing experiments, SEQC.py provides some convenience tools to create indices for use with SEQC and STAR, and tools to check the progress of remote runs, list current runs, start instances, and terminate them. 
