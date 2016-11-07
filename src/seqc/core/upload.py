import os
import pandas as pd
from seqc import io, log
import shutil
from subprocess import Popen, PIPE
from seqc.stats.experimental_yield import ExperimentalYield
from seqc.ec2 import Retry


# todo make this function less complex
def data_and_notify(
        email_status: str, aws_upload_key: str, align: bool, input_data: str,
        manage_merged: io.ProcessManager, process_samfile: bool,
        merged_fastq: str, samfile: str, read_array: str,
        sparse_proc: io.ProcessManager, sparse_csv: str,
        manage_ra: io.ProcessManager, summary: dict, fastq_records: int,
        sam_records: int, total_molecules: int, mols_lost: int, cells_lost: int,
        manage_samfile: io.ProcessManager, cell_description: pd.Series,
        output_stem: str, log_name: str, email: bool) -> None:
    """
    Uploads data and notifies the user of the termination of a run

    :param email_status: str, email address to email results
    :param aws_upload_key: str, location to upload files on aws
    :param align: bool, whether or not alignment occurred
    :param input_data: str, indicator, stores the type of input data
    :param manage_merged: io.ProcessManager, an upload manager for merged_fastq
    :param process_samfile: bool, whether or not samfile was processed
    :param merged_fastq: str, name of merged fastq file
    :param samfile: str, name of sam file
    :param read_array: str, name of stored ReadArray
    :param sparse_proc: io.ProcessManager, an upload manager for sparse counts file
    :param sparse_csv: str, name of sparse_csv file
    :param manage_ra: io.ProcessManager, an upload manager for the ReadArray
    :param summary: dict, object containing summary statistics
    :param fastq_records: int, number of processed fastq records
    :param sam_records: int, number of processed sam records
    :param total_molecules: int, number of molecules
    :param mols_lost: int, number of molecules lost
    :param cells_lost: int, number of cells lost
    :param manage_samfile: io.ProcessManager, an upload manager for the sam file
    :param cell_description: pd.Series, a summary of cells and cell counts
    :param output_stem: str, stem for file output
    :param log_name: str, name of SEQC.log file
    :param email: bool, whether email should be sent
    :return None:
    """
    
    #TODO: Some parts of it are commented out as they rely on stuff returned from "generate_count_matrices()" which is now deprecated and replaced.
    # Need to think what is missing from the summary now and add it.

    log.info('Starting file upload onto %s.' % aws_upload_key)

    if email_status and aws_upload_key is not None:
        # make sure that all other files are uploaded before termination
        if align and input_data != 'merged':
            if manage_merged:
                manage_merged.wait_until_complete()
                log.info('Successfully uploaded %s to the specified S3 location "%s"' %
                         (merged_fastq, aws_upload_key))
        if process_samfile:
            if input_data != 'samfile':
                if manage_samfile:
                    manage_samfile.wait_until_complete()
                    log.info('Successfully uploaded %s to the specified S3 location "%s"'
                             % (samfile, aws_upload_key))
            if manage_ra:
                manage_ra.wait_until_complete()
                log.info('Successfully uploaded %s to the specified S3 location "%s"' %
                         (read_array, aws_upload_key))
            #sparse_proc.wait_until_complete()
            #log.info('Successfully uploaded %s to the specified S3 location "%s"' %
            #         (sparse_csv + '.gz', aws_upload_key))

        # upload count matrix and alignment summary at the very end
        if summary:
            if fastq_records:
                summary['n_fastq'] = fastq_records
            else:
                summary['n_fastq'] = 'NA'
            if sam_records:
                summary['n_sam'] = sam_records
            else:
                summary['n_sam'] = 'NA'
            #summary['total_mc'] = total_molecules
            #summary['mols_lost'] = mols_lost
            #summary['cells_lost'] = cells_lost
            #summary['cell_desc'] = cell_description

        # todo: run summary will not be reported if n_fastq or n_sam = NA
        upload_results(
            output_stem, email_status, aws_upload_key, input_data,
            summary, log_name, email)


def upload_results(output_stem: str, email_address: str, aws_upload_key: str,
                   start_pos: str, summary: dict, log_name: str, email: bool) -> None:
    """
    uploads remaining files from the SEQC run
    :param output_stem: specified output directory in cluster
    :param email_address: e-mail where run summary will be sent
    :param aws_upload_key: tar gzipped files will be uploaded to this S3 bucket
    :param start_pos: determines where in the script SEQC started
    :param summary: dictionary of summary statistics from SEQC run
    :param log_name: log name of SEQC run provided by user
    :param email: True if mutt is installed, False otherwise.
    """

    prefix, directory = os.path.split(output_stem)
    if not prefix:
        prefix = '.'
    counts = output_stem + '_read_and_count_matrices.p'
    log_text = prefix + '/' + log_name

    # generate a run summary and append to log + email
    run_summary = ExperimentalYield.construct_run_summary(summary)
    log.info('A copy of the SEQC run summary can be found below.\nRUN SUMMARY:\n'
             '{run_summary}'.format(run_summary=run_summary))
    files = [counts, log_text]  # counts and log will always be uploaded

    if start_pos == 'start' or start_pos == 'merged':
        alignment_summary = output_stem + '_alignment_summary.txt'
        # copying over alignment summary for upload
        shutil.copyfile(prefix + '/alignments/Log.final.out', output_stem +
                        '_alignment_summary.txt')
        files.append(alignment_summary)

    if aws_upload_key:
        bucket, key = io.S3.split_link(aws_upload_key)
        for item in files:
            try:
                Retry(retries=5)(io.S3.upload_file)(item, bucket, key)
                item_name = item.split('/')[-1]
                log.info('Successfully uploaded %s to the specified S3 location '
                         '"%s%s".' % (item, aws_upload_key, item_name))
            except FileNotFoundError:
                log.notify('Item %s was not found! Continuing with upload...' % item)

        # get the name of the output file
        log.info('Upload complete. An e-mail will be sent to %s.' % email_address)

    # email results to user
    body = ('<font face="Courier New, Courier, monospace">'
            'SEQC RUN COMPLETE.\n\n'
            'The run log has been attached to this email and '
            'results are now available in the S3 location you specified: '
            '"%s"\n\n'
            'RUN SUMMARY:\n\n%s'
            '</font>' % (aws_upload_key, run_summary))
    log.info('SEQC run complete. Cluster will be terminated unless --no-terminate '
             'flag was specified.')
    body = body.replace('\n', '<br>')
    body = body.replace('\t', '&emsp;')
    if email:
        email_user(log_text, body, email_address)


# todo refactor, needs its own module
def email_user(attachment: str, email_body: str, email_address: str) -> None:
    """
    sends an email to email address with text contents of email_body and attachment
    attached. Email will come from "Ubuntu@<ec2-instance-ip-of-aws-instance>

    :param attachment: the file location of the attachment to append to the email
    :param email_body: text to send in the body of the email
    :param email_address: the address to which the email should be sent"""

    if isinstance(email_body, str):
        email_body = email_body.encode()

    email_args = ['mutt', '-e', 'set content_type="text/html"', '-a', attachment, '-s',
                  'Remote Process', '--', email_address]
    email_process = Popen(email_args, stdin=PIPE)
    email_process.communicate(email_body)
