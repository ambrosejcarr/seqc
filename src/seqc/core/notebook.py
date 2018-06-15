from seqc.notebooks.notebooks import Notebook
from seqc import log


def notebook(args):
    if args.subsubparser_name == 'merge':
        # need to also take a output directory because this thing will write stuff.
        # then merge the things
        # then return?
        n = Notebook(args.output_filename, *args.input_data)
        n.merge_data(merged_sample_name=args.output_filename)
        log.info('Merged samples written to %s' % args.input_data)
    elif args.subsubparser_name == 'generate':
        n = Notebook(args.output_stem, args.input_count_matrix)
        n.write_template()
        log.info('Notebook Template written to %s' % n.notebook_path)
        n.run_notebook()
        log.info('Notebook Run and written to %s' % n.notebook_path)

