import gzip


def create_gtf(annotation_file, nums_to_transcripts_file, gtf_out_file):
    if annotation_file.endswith('.gz'):
        in_file = gzip.open(annotation_file, 'rt')
    else:
        in_file = open(annotation_file)
    in_read = in_file.readlines()
    in_file.close()

    nums_file = open(nums_to_transcripts_file)
    nums_read = nums_file.readlines()
    nums_file.close()

    # So the goal is just to turn the annotations.gtf into a new file that will be on our
    # scseq_id feature level. We have the following fields to consider:
    # [CDS, exon, gene, transcript, UTR] and possibly also [Selenocysteine,
    # start_codon, stop_codon]
    # Probably the last two in the second category but maybe not Selenocysteine unless
    # necessary? Only 60 sites for it across the whole thing...

    tran_to_num = {}
    num_to_tran = {}
    for line in nums_read:
        lsp = line.split(",")
        num_to_tran[int(lsp[0].split("SC")[1])] = lsp[1].strip()
        tran_to_num[lsp[1].strip()] = int(lsp[0].split("SC")[1])

    out_file = open(gtf_out_file, "w")
    for line in in_read:
        if line[0] == "#":
            out_file.write(line)
        else:
            lsp = line.split("\t")
            if lsp[2] != "gene" and "transcript_id" in lsp[-1]:
                tran_id = \
                    lsp[-1].split("transcript_id")[1].split(";")[0].split('"')[1].split('.')[
                        0]
                gene_id = \
                    lsp[-1].split("gene_id")[1].split(";")[0].split('"')[1].split('.')[0]
                try:
                    key = tran_id + "_" + gene_id
                    out_file.write(
                        line.strip() + ' scseq_id "SC' + str(tran_to_num[key]) + '";\n')
                except KeyError:
                    #print tran_id + "_" + gene_id
                    out_file.write(line.strip() + ' scseq_id "SC000";\n')
                #missing.append(tran_id + "_" + gene_id)
            else:
                out_file.write(line)
    out_file.close()

# err_file.write("\n".join(missing) + '\n')
#err_file.close()

if __name__ == "__main__":
    create_gtf("/root/home/sa_preprocess/annotations.gtf",
               "/mnt/jelly/test_nums_to_transcripts.txt", "reduced_gtf_w_gene.gtf")
