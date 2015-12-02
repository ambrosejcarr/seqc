#!/usr/local/bin/python3

from subprocess import Popen, PIPE, check_output
import seqc
import argparse


def main(cmdstring):
#cmdstring = "SEQC in-drop --forward short_f1.fastq --index s3://dplab-data/genomes/mm38/ --frag-len 1000 --reverse short_r1.fastq --n-threads 30 --barcodes s3://dplab-data/sc-seq/allon/barcodes/in_drop_barcodes.p --output-file /data/software/qq"
# cmdstring = "SEQC in-drop --forward fx.fastq --index s3://dplab-data/genomes/mm38/ --frag-len 1000 --reverse short_r1.fastq --n-threads 30 --barcodes s3://dplab-data/sc-seq/allon/barcodes/in_drop_barcodes.p --output-file /data/software/qq"

    cmd = cmdstring.split(' ')
    p = Popen(cmd, stderr=PIPE, stdout=PIPE)
    out, err = p.communicate()

    if err:
        errmsg = err.decode().split('\n')
        with open('errors.txt', 'w') as f:
            for x in errmsg:
                f.write(x + '\n')
        ps = Popen(['echo', "Process interrupted -- see attached error message"],
                   stdout=PIPE, stderr=PIPE)
        output = Popen(['mutt', '-a', 'errors.txt', '-s', 'Remote Process', '--',
                        'kristy.choi24@gmail.com'],
                        stdin=ps.stdout, stdout=PIPE, stderr=PIPE)
    else:
        ps = Popen(['find', '.', '-type', 'f'], stdout=PIPE, stderr=PIPE)
        output = Popen(['grep', '-E', 'out.sam$|.h5$|merged|.npz$'], stdin=ps.stdout,
                       stdout=PIPE, stderr=PIPE)

        # get name of .npz file here
        pp = Popen(['find', '.', '-type', 'f'], stdout=PIPE, stderr=PIPE)
        pp_out = Popen(['grep', '.npz$'], stdin=pp.stdout, stdout=PIPE, stderr=PIPE)
        out, err = pp_out.communicate()
        res = out.decode().split('\n')[0]

        # make upload directory and gzip resulting files
        makedir = check_output(['mkdir', 'upload'])
        fin = check_output(['xargs', 'cp', '-t', 'upload'], stdin=output.stdout)
        zipped = check_output(['tar', '-czf', 'seqc_out.tar.gz', 'upload'])
        seqc.io.S3.upload_file('seqc_out.tar.gz', 'dplab-data', 'seqc/test_seqc/')

        # in the case that .npz file exists, which should be always
        if res:
            ms = Popen(['echo',
                    "seqc run complete -- see attached .npz file. The rest of the output is "
                    "available as seqc_out.tar.gz in "
                    "your specified S3 bucket."], stdout=PIPE, stderr=PIPE)
            msgout = Popen(['mutt', '-a', res, '-s', 'Remote Process', '--',
                            'kristy.choi24@gmail.com'], stdin=ms.stdout,
                       stdout=PIPE, stderr=PIPE)
        else:
            # just for now when we don't have the .npz file
            ms = Popen(['echo', "seqc run complete! The rest of the output is available as "
                                "seqc_out.tar.gz in "
                            "your specified S3 bucket."], stdout=PIPE, stderr=PIPE)
            msgout = Popen(['mutt', '-s', 'Remote Process', '--', 'kristy.choi24@gmail.com'],
                           stdin=ms.stdout, stdout=PIPE,
                           stderr=PIPE)

if __name__ == "__main__":
    #process cmdstring
    p = argparse.ArgumentParser()
    p.add_argument('-l', help = 'SEQC call')
    args = p.parse_args()

    #run script
    main(args.l)
