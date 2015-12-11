

# of all the commands in core, only merged_fastq takes subparser_name
# merge converts that into exp_type and uses it for the 3bitprocessor

# other methods that use the processors are:
# (1) Tests
# (2) SEQC script -- needs to know about barcodes (y/n) and earliest insertion into
#     pipeline; currently has modules for each experiment type, but we could change that
#     to factory classes/functions
# (3) barcodes.CellBarcodes is aware of DropSeq, but others it has a standard form.
# (4) encodings.ThreeBit
# (5) generator functions for fastq
# (6) generator functions for sam


# idea: the below class, or .ini file or function or dictionary or however we want to
# define it gets called by the SEQC script to define the runtime environment. Everything
# needed to create a run should therefore be present in the class. For things that are
# just too complicated to get wrapped up, we can define custom pre-processing methods
# and add them to be pre-run.


class DropSeq():

    # which arguments to add to parser? if not forward/reverse, then merge will not
    #  be run. If not merged, then align will not be run.. etc. Builds functions to
    #  process the data dynamically.
    # Note that below 'forward' and 'reverse' have been replaced with 'genomic' and
    #  'barcodes', to indicate which strand contains barcode data (some are reversed).
    valid_inputs = ['genomic', 'barcodes', 'merged', 'sam']

    # if barcodes_required is false, evidence-based error correction will not be run
    barcodes_required = False

    # methods to be run prior to beginning of general pipeline (which runs from merged
    # fastq); if present, this method replaces the general CellBarcodes method, which
    # should take a forward and reverse read and merge them into the standard form.
    custom_merge_method = None

    # general preprocessing parameters; these data should be adequate to create a
    # generator method for each data type.
    cell_barcode_position = (0, 12)
    umi_position = (12, 20)



