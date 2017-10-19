import pandas as pd
import numpy as np
import gzip
import os
from seqc.sequence.fastq import Reader
from seqc.sequence.barcodes import hamming_dist_bin
from seqc.sequence.encodings import DNA3Bit
from seqc import log, ec2, platforms, io
from copy import deepcopy
import time
from seqc.core import verify, download
from seqc.email_ import email_user
import re
import gzip
    
def demultiplex(args):
    
    def _read_sample_list_file(sample_index_file):
        sample_indexes = []
        sample2outputprefix = dict()
        file = open(sample_index_file, 'r')
        for line in file.readlines():
            arr = re.split("[\t ,]+", line)
            arr[-1] = arr[-1].rstrip()
            sample_indexes.append(arr[0])
            sample2outputprefix[arr[0]] = arr[1]
        file.close()
        return sample_indexes, sample2outputprefix
    
    def _read_index_file(index_file):
        indexes = pd.read_csv(index_file, header=None, index_col=0)
    
        # Mapping from indexes to samples 
        indexes = pd.Series(np.repeat(indexes.index, indexes.shape[1]),
            index = [DNA3Bit.encode(i.encode()) for i in np.ravel(indexes)])
    
        return indexes
    
    def _write_record(r, handle ):
        # Write record
        handle.writelines([r.name, r.sequence, r.name2, r.quality])


    def demultiplexSamples(fq_files, index_file, sample_indexes,
                    output_dir, sample2outputprefix, max_ed=1):

        output_files=[]                 #    for returning the list of output files
        outputfile2sample=dict()        #    mapping each output file to sample index 
        
        # Read indexes 
        valid_indexes = _read_index_file(index_file)
    
        # Restrict to valid indexes
        valid_indexes = valid_indexes[valid_indexes.isin(sample_indexes)]
        index_list = list(valid_indexes.index)
        # Convert to dictionary for speedup
        valid_indexes = valid_indexes.to_dict()
    
    
        # Counters
        counters = dict()
        counters['CORRECT_INDEX'] = 0
        counters['CORRECTED_INDEX'] = 0
        counters['INCORRECT_INDEX'] = 0
        sample_counters = dict()
        
        sample2outputprefix['Undetermined']="Sample"
    
        # Open file handles
        file_handles = dict()
        for sample in np.append(list(set(valid_indexes.values())), ['Undetermined']):
            sample_counters[sample] = deepcopy(counters)
            for r in ['R1', 'R2']:
                filename = output_dir + "/" + sample2outputprefix[sample] + "_" + sample + '_' + r + '.fastq'+".gz"
                output_files.append(output_dir + "/" + sample2outputprefix[sample] + "_" +sample + '_' + r + '.fastq'+".gz")
                outputfile2sample[output_dir + "/" + sample2outputprefix[sample] + "_" +  sample + '_' + r + '.fastq'+".gz"] = sample
                file_handles[filename] = gzip.open(filename, 'wb')
    
        # PolyN counter
        polyn = DNA3Bit.encode(b'NNNNNNNN')
        sample_counters[polyn] = 0
    
        sampled=dict()    #    save already computed sample for already seen index
        mdistd=dict()     #    save already computed minimum distance for already seen index  
    
        # Forward and reverse iterators
        numr = 0
        for f, r in zip(Reader(fq_files[0]), Reader(fq_files[1])):
    
            # Extract index
            index = DNA3Bit.encode(f.name.strip().split(b':')[::-1][0].split(b'+')[0])
    
            # Skip past polyN
            if index == polyn:
                sample_counters[polyn] += 1
                continue
    
            # Hamming distance to all indexes
            if index in sampled:
                sample = sampled[index]
                mdist = mdistd[index]
            else:
                dists = [hamming_dist_bin(index, i) for i in index_list]
                mdist = np.min(dists)
                if mdist > max_ed:
                    sample = 'Undetermined'
                else:
                    sample = valid_indexes[index_list[np.argmin(dists)]]
                
                mdistd[index] = mdist
                sampled[index] = sample
                            
            # Counting the indexes
            if sample == 'Undetermined':
                sample_counters[sample]['INCORRECT_INDEX'] += 1
            else:
                if mdist == 0:
                    sample_counters[sample]['CORRECT_INDEX'] += 1
                else:
                    sample_counters[sample]['CORRECTED_INDEX'] += 1
            # Write tofile
            
            
            
            _write_record(f, file_handles[output_dir + "/" + sample2outputprefix[sample] + "_" + sample + '_R1.fastq'+".gz"])
            _write_record(r, file_handles[output_dir + "/" + sample2outputprefix[sample] + "_" + sample + '_R2.fastq'+".gz"])
            
            numr = numr + 1
            #if numr>1000:
                #break
            
        # Close the file handles
        for k in file_handles.keys():
            file_handles[k].close()
    
        return sample_counters, output_files, outputfile2sample
    
    log.setup_logger(args.log_name)

    with ec2.instance_clean_up(
            email=args.email, upload=args.upload_prefix, log_name=args.log_name,
            debug=args.debug):
        pigz, mutt = verify.executables('pigz', 'mutt')
        if mutt:
            log.notify('mutt executable identified, email will be sent when run '
                       'terminates. ')
        else:
            log.notify('mutt was not found on this machine; an email will not be sent to '
                       'the user upon termination of SEQC run.')
        log.args(args)
        
        output_dir = '.'
            
        # download index-map-file if has to
        if args.index_map_file.startswith("s3:/"):
            download.s3_data(
                [args.index_map_file],output_dir)
            index_map_file=output_dir+"/"+os.path.basename(args.index_map_file)
        
        # download sample-list-file if has to    
        if args.sample_list_file.startswith("s3:/"):
            download.s3_data([args.sample_list_file],output_dir)
            sample_list_file=output_dir+"/"+os.path.basename(args.sample_list_file)
        
        log.info("downloaded sample list file and index map file")
     
        # downloading fastq files
        log.info("downloading pooled fastq files")
        genomic_fastq = download.s3_data(args.genomic_fastq, output_dir+"/genomic_fastq/")
        barcode_fastq = download.s3_data(args.barcode_fastq, output_dir+"/barcode_fastq/")
        
        sample_indexes, sample2outputprefix=_read_sample_list_file(sample_list_file)
        log.info("read sample list file")
        
        output_files = []
        log.info("starting demultiplexing")
        sample_counters, output_files, outputfile2sample=demultiplexSamples([barcode_fastq,genomic_fastq], index_map_file, sample_indexes, output_dir, sample2outputprefix)
        log.info("demultiplexing finished")
        
        polyn = DNA3Bit.encode(b'NNNNNNNN')
        log.info("SUMMARY REPORT:")
        log.info("%25s\t%25s\t%25s" % ("SAMPLE_INDEX", "#EXACT", "#CORRECTED"))
        for s in sample2outputprefix:
            if (s != 'Undetermined') and (s != polyn):
                log.info("%25s\t%25s\t%25s" % (s, str(sample_counters[s]['CORRECT_INDEX']), str(sample_counters[s]['CORRECTED_INDEX'])))
        log.info("%25s\t%25s" % ("Undetermined", str(sample_counters['Undetermined']['INCORRECT_INDEX'])))
        log.info("%25s\t%25s" % ('NNNNNNNN', str(sample_counters[polyn])))
        
        upload_files = output_files
        upload_files.append(args.log_name)
        upload_files.append('./nohup.log')
        
        
        if args.upload_prefix:
            
            # Upload the sample fastq files and log files            
            for item in upload_files:
                if "R1.fastq.gz" in item:
                    sampleIndex = outputfile2sample[item]
                    bucket, key = io.S3.split_link(args.upload_prefix + sample2outputprefix[sampleIndex]+ "_"  + sampleIndex + "/barcode_fastq/")
                elif "R2.fastq.gz" in item:
                    sampleIndex = outputfile2sample[item]
                    bucket, key = io.S3.split_link(args.upload_prefix + sample2outputprefix[sampleIndex]+ "_"  + sampleIndex + "/genomic_fastq/")
                else:
                    bucket, key = io.S3.split_link(args.upload_prefix)
                    
                try:
                    ec2.Retry(retries=5)(io.S3.upload_file)(item, bucket, key)
                    item_name = item.split('/')[-1]
                    log.info('Successfully uploaded %s to the specified S3 location '
                             '"%s%s".' % (item, "s3://"+bucket+"/"+key, item_name))
                except FileNotFoundError:
                    log.notify('Item %s was not found! Continuing with upload...' % item)

        # todo local test does not send this email
        if mutt:
            email_body = (
                '<font face="Courier New, Courier, monospace">'
                'SEQC RUN COMPLETE.\n\n'
                'The run log has been attached to this email and '
                'results are now available in the S3 location you specified: '
                '"%s"\n\n' % args.upload_prefix)
            email_body = email_body.replace('\n', '<br>').replace('\t', '&emsp;')
            email_user("", email_body, args.email)