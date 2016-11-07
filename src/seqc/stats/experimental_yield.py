class ExperimentalYield:

    output = (
        '{divide}\nINPUT\n{divide}\n'
        'Total input reads:\t{n_fastq}\n'
        '{divide}\nALIGNMENT (% FROM INPUT)\n{divide}\n'
        'Total reads aligned:\t{n_sam} ({prop_al}%)\n'
        ' - Genomic alignments:\t{genomic} ({prop_gen}%)\n'
        ' - PhiX alignments:\t{phi_x} ({prop_phix}%)\n'
        ' - Transcriptome alignments:\t{trans} ({prop_trans}%)\n'
        '{divide}\nFILTERING (% FROM ALIGNMENT)\n{divide}\n'
        'Genomic alignments:\t{genomic} ({bad_gen}%)\n'
        'PhiX alignments:\t{phi_x} ({bad_phi}%)\n'
        'Incorrect barcodes:\t{wrong_cb} ({bad_cb}%)\n'
        'Missing cell barcodes/RMT:\t{no_cell} ({bad_cell}%)\n'
        'N present in RMT:\t{rmt_N} ({bad_rmtN}%)\n'
        'N present in CB:\t{cell_N} ({bad_cellN}%)\n'
        'Insufficient poly(T):\t{poly_t} ({bad_polyt}%)\n'
        'High dust score:\t{dust} ({bad_dust}%)\n'
        '{divide}\nCELL/MOLECULE COUNT DISTRIBUTION\n{divide}\n'
        'Total molecules:\t\t{tot_mc}\n'
        'Molecules lost:\t{mols_lost}\n'
        'Cells lost:\t{cells_lost}\n'
        'Cell description:\n{cell_desc}\n'
        '{divide}\nSUMMARY\n{divide}\n'
        'Total retained reads:\t{n_good} ({prop_good}%)\n'
        'Total reads unaligned:\t{lost_al} ({prop_un}%)\n'
        'Total reads filtered:\t{n_bad} ({prop_bad}%)\n'
        '{divide}\n')

    @classmethod
    def construct_run_summary(cls, summary: dict):
        """
        calculates basic loss statistics and constructs a summary
        that will be sent to the user after the SEQC run has completed.

        :param summary: dictionary constructed during error correction
        :return: output of basic summary statistics
        """
        if not summary:
            return

        # obtain values from summary
        n_fastq = summary['n_fastq']
        n_sam = summary['n_sam']
        genomic = summary['gene_0']
        phix = summary['phi_x']
        no_cell = summary['cell_0']
        # no_rmt = summary['rmt_0']
        rmt_N = summary['rmt_N']
        cell_N = summary['cell_N']
        dust = summary['dust']
        poly_t = summary['poly_t']
        tot_mc = summary['total_mc']
        mols_lost = list(summary['mols_lost'].items())
        cells_lost = list(summary['cells_lost'].items())
        cell_desc = summary['cell_desc'].to_string()
        divide = '-' * 40

        # run summary will not be calculated if user started SEQC midway
        if n_fastq == 'NA' or n_sam == 'NA':
            return

        # calculate summary statistics
        trans = n_sam - genomic - phix
        prop_al = round((n_sam/n_fastq) * 100, 1)
        prop_gen = round((genomic/n_sam) * 100, 1)
        prop_phix = round((phix/n_sam) * 100, 1)
        prop_trans = round((trans/n_sam) * 100, 1)
        lost_al = n_fastq - n_sam
        prop_un = round(100 - prop_al, 1)
        n_bad = genomic + phix + no_cell + rmt_N + cell_N + poly_t + dust
        # n_bad = genomic + phix + no_cell + no_rmt + rmt_N + poly_t
        # wrong_cb does not apply to drop-seq
        try:
            wrong_cb = summary['cb_wrong']
            n_bad += wrong_cb
            bad_cb = round((wrong_cb/n_bad) * 100, 1)
        except KeyError:
            wrong_cb = 0
            bad_cb = 0
        # continue with calculations
        n_good = n_sam - n_bad
        bad_gen = round((genomic/n_bad) * 100, 1)
        bad_phi = round((phix/n_bad) * 100, 1)
        bad_cell = round((no_cell/n_bad) * 100, 1)
        # bad_rmt = round((no_rmt/n_bad) * 100, 1)
        bad_rmtN = round((rmt_N/n_bad) * 100, 1)
        bad_cellN = round((cell_N/n_bad) * 100, 1)
        bad_polyt = round((poly_t/n_bad) * 100, 1)
        bad_dust = round((dust/n_bad) * 100, 1)
        prop_bad = round((n_bad/n_fastq) * 100, 1)
        prop_good = round((n_good/n_fastq) * 100, 1)

        # format output
        output = cls.output.format(
            n_fastq=n_fastq, n_sam=n_sam, genomic=genomic, phi_x=phix, no_cell=no_cell,
            wrong_cb=wrong_cb, rmt_N=rmt_N, poly_t=poly_t, divide=divide,
            prop_al=prop_al, prop_gen=prop_gen, prop_phix=prop_phix, lost_al=lost_al,
            n_bad=n_bad, n_good=n_good, prop_good=prop_good, prop_bad=prop_bad,
            prop_un=prop_un, bad_gen=bad_gen, bad_phi=bad_phi, bad_cb=bad_cb,
            bad_cell=bad_cell, bad_rmtN=bad_rmtN, bad_polyt=bad_polyt, trans=trans,
            cell_N=cell_N, bad_cellN=bad_cellN, dust=dust, bad_dust=bad_dust,
            prop_trans=prop_trans, tot_mc=tot_mc, mols_lost=mols_lost,
            cells_lost=cells_lost, cell_desc=cell_desc)
        return output
