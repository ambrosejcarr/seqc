# A single file, called prepare_fasta, is needed to produce some files that will be used
# in downstream steps of suffix array construction.

from seqc.sa_preprocess import load_fasta as lf
from seqc.sa_preprocess import buffered_printer as bf
from seqc.log import log_info
import random
import gzip
import pickle


def random_nuc():
    """Returns random nucleotide (for remove_Ns)."""
    return 'ACGT'[int(random.random() * 4)]


def remove_Ns(seq):
    """Replaces Ns with random nucleotides. Ns wouldn't be considered potential overlaps
    and needlessly complicate suffix array, and are present in only 6 spots in mouse
    transcriptome. Prefer to have very small risk of false-positives from randomly adding
    a letter.
    """
    while "N" in seq:
        pos_to_change = seq.find("N")
        if pos_to_change < len(seq) - 1:
            seq = seq[:pos_to_change] + random_nuc() + seq[pos_to_change + 1:]
        else:
            seq = seq[:pos_to_change] + random_nuc()
    return seq


# This code is needed as it is used by the code in the ensembl_76 directory,
# new_filter_ensembl_cdna_and_ncrna.py

def find_acceptable(annotations_gtf, fasta_file):
    """Acceptable transcripts are described in the genome's annotation file (gtf) and in
    the genome's fasta file, so that we can descriptions of genomic context as well as
    sequence when combining those two sources. We need both. This returns list of ids
    for 'acceptable' transcripts (ids will be the header line for each sequence in the
    fasta file)."""

    if annotations_gtf.endswith('.gz'):
        in_file = gzip.open(annotations_gtf, 'rt')
    else:
        in_file = open(annotations_gtf)
    in_read = in_file.readlines()
    in_file.close()

    trans = []
    for line in in_read:
        if line.startswith("#"):
            continue
        else:
            lsp = line.split("\t")
            if lsp[2] != "transcript":
                continue
            meta = lsp[8]
            tran_id = [field for field in meta.split(";") if "transcript_id" in field]
            gene_id = [field for field in meta.split(";") if "gene_id" in field]
            assert len(tran_id) <= 1
            assert len(gene_id) <= 1
            tran_id = tran_id[0].split("transcript_id ")[1].strip('"').split(".")[0]
            gene_id = gene_id[0].split("gene_id ")[1].strip('"').split(".")[0]
            trans.append(tran_id)

    trans = set(trans)  # the set of all transcript ids...

    if fasta_file.endswith('.gz'):
        my_file = gzip.open(fasta_file, 'rt')
    else:
        my_file = open(fasta_file)
    my_read = my_file.readlines()
    my_file.close()

    my_dict = lf.load_fasta(my_read)
    my_keys = list(my_dict.keys())

    accepted_dict = {}
    for k in my_dict:
        if k.split()[0] in trans:
            accepted_dict[k] = my_dict[k]
    return accepted_dict


def process_sequences(accepted_dict, out_file_name):
    """
    Sequences have to be processed so that they don't contain Ns and so that they are
    the correct length."""

    acceptable_transcripts = {}
    for k in accepted_dict:
        cur_gene = accepted_dict[k]
        if "N" in cur_gene.upper():
            cur_gene = remove_Ns(cur_gene.upper())
        acceptable_transcripts[k] = cur_gene

    accepted_dict = {}

    to_output = {}
    gene_termini = {}

    at_counter = 0

    all_tails = {}
    all_filter = {}

    for at in acceptable_transcripts:
        at_counter += 1
        tail = acceptable_transcripts[at][-1000:]
        gene = at.split(" gene:")[1].split(" ")[0]
        tran = at.split(" ")[0]
        to_output[tran + "_" + gene] = tail

    out_file = open(out_file_name, "w")
    for tran_id in to_output:
        out_file.write(
            ">" + tran_id + "\n" + bf.buffered_printer(to_output[tran_id], 60) + "\n")
    out_file.close()


def standard_run(annotations_gtf="annotations.gtf",
                 cdna_plus_ncrna_fa="Mus_musculus.GRCm38.cdna_all_plus_ncrna.fa",
                 output_file="csf_labeled_fasta.fa"):
    """Takes arguments: annotations_gtf, cdna_plus_ncrna_fa, output_file. All paths to
    file (not file objects). The file annotations_gtf is the gtf file giving gene
    annotations for the organism of interest. cdna_plus_ncrna_fa is the concatenated file
    of cdna.fa and ncrna.fa for organism of interest. output_file is the desired local
    path of the file to be outputted by this program.
    """
    log_info("Preprocessing fasta and gtf files")
    accepted_dict = find_acceptable(annotations_gtf, cdna_plus_ncrna_fa)
    process_sequences(accepted_dict, output_file)


if __name__ == "__main__":
    standard_run("annotations.gtf", "Mus_musculus.GRCm38.cdna_all_plus_ncrna.fa",
                 "csf_labeled_fasta.fa")
