import regex as re


_pattern = re.compile(b'(.{8,11}?)(GAGTGATTGCTTGTGACGCCTT){s<=2}(.{8})(.{6})(.*?)')
_pattern_v2 = re.compile(b'(.{8,11}?)(GAGTGATTGCTTGTGACGCCAA){s<=2}(.{8})(.{8})(.*?)')


def _check_spacer_v2(sequence):
    """a fast, in-drop-v2 specific command to find a spacer sequence and cb length

    :param sequence: fastq sequence data
    :returns: (cell_barcode, rmt, poly_t)
    """
    assert ~sequence.endswith(b'\n')
    identifier = sequence[24:28]
    if identifier == b'CGCC':
        cb1 = sequence[:8]
        cb2 = sequence[30:38]
        rmt = sequence[38:46]
        poly_t = sequence[46:]
    elif identifier == b'ACGC':
        cb1 = sequence[:9]
        cb2 = sequence[31:39]
        rmt = sequence[39:47]
        poly_t = sequence[47:]
    elif identifier == b'GACG':
        cb1 = sequence[:10]
        cb2 = sequence[32:40]
        rmt = sequence[40:48]
        poly_t = sequence[48:]
    elif identifier == b'TGAC':
        cb1 = sequence[:11]
        cb2 = sequence[33:41]
        rmt = sequence[41:49]
        poly_t = sequence[49:]
    else:
        return b'', b'', b''
    cell = cb1 + cb2
    return cell, rmt, poly_t


def _check_spacer_v1(sequence):
    """a fast, in-drop-v1 specific command to find a spacer sequence and cb length

    :param sequence: fastq sequence data
    :returns: (cell_barcode, rmt, poly_t)
    """
    assert ~sequence.endswith(b'\n')
    identifier = sequence[24:28]
    if identifier == b'CGCC':
        cb1 = sequence[:8]
        cb2 = sequence[30:38]
        rmt = sequence[38:44]
        poly_t = sequence[44:]
    elif identifier == b'ACGC':
        cb1 = sequence[:9]
        cb2 = sequence[31:39]
        rmt = sequence[39:45]
        poly_t = sequence[45:]
    elif identifier == b'GACG':
        cb1 = sequence[:10]
        cb2 = sequence[32:40]
        rmt = sequence[40:46]
        poly_t = sequence[46:]
    elif identifier == b'TGAC':
        cb1 = sequence[:11]
        cb2 = sequence[33:41]
        rmt = sequence[41:47]
        poly_t = sequence[47:]
    else:
        return b'', b'', b''
    cell = cb1 + cb2
    return cell, rmt, poly_t


def in_drop(g, b):
    """
    merge forward and reverse in-drop v1 reads, annotating the reverse read (containing
    genomic information) with the cell barcode, rmt, and number of poly_t. Pool is left
    empty.

    :param g: genomic fastq sequence data
    :param b: barcode fastq sequence data
    :return: annotated genomic sequence.
    """

    cell, rmt, poly_t = _check_spacer_v1(b.sequence[:-1])
    if not cell:
        try:
            cell1, spacer, cell2, rmt, poly_t = re.match(
                _pattern, b.sequence[:-1]).groups()
            cell = cell1 + cell2
        except AttributeError:
            cell, rmt, poly_t = b'', b'', b''
    g.add_annotation((b'', cell, rmt, poly_t))
    return g


def in_drop_v2(g, b):
    """
    merge forward and reverse in-drop v2 reads, annotating the reverse read (containing
    genomic information) with the cell barcode, rmt, and number of poly_t. Pool is left
    empty.

    :param g: genomic fastq sequence data
    :param b: barcode fastq sequence data
    :return: annotated genomic sequence.
    """
    cell, rmt, poly_t = _check_spacer_v2(b.sequence[:-1])
    if not cell:
        try:
            cell1, spacer, cell2, rmt, poly_t = re.match(
                _pattern_v2, b.sequence[:-1]).groups()
            cell = cell1 + cell2
        except AttributeError:
            cell, rmt, poly_t = b'', b'', b''
    g.add_annotation((b'', cell, rmt, poly_t))
    return g


def in_drop_v3(g, b):
    """
    merge forward and reverse in-drop v3 reads, annotating the reverse read (containing
    genomic information) with the rmt and number of poly_t from the
    forward read. Pool is left empty, and the cell barcode is reconstructed from the
    second index and the second barcode.

    Please note that R1 is genomic, and R2 is barcode, unlike previous iterations

    :param g: genomic fastq sequence data
    :param b: barcode fastq sequence data
    :return: annotated genomic sequence.
    """
    seq = b.sequence.strip()
    cell2 = seq[:8]
    rmt = seq[8:16]
    poly_t = seq[16:]
    # bc is in a fixed position in the name; assumes 8bp indices.
    cell1 = g.name.strip()[-17:-9]
    g.add_annotation((b'', cell1 + cell2, rmt, poly_t))
    return g


def ten_x(g, b):
    """
    merge forward and reverse 10x reads, annotating the reverse read
    (containing genomic information) with the rmt from the forward read.
    Pool is left empty, and the cell barcode is obtained from the
    name field of the forward read.

    Please note that R1 is genomic, and R2 is RMT, unlike previous iterations

    :param g: genomic fastq sequence data
    :param b: barcode fastq sequence data
    :return: annotated genomic sequence.
    """
    rmt = b.sequence.strip()
    # bc is in a fixed position in the name; assumes 10bp indices.
    cell = g.name.strip()[-23:-9]
    g.add_annotation((b'', cell, rmt, b''))
    return g


def drop_seq(g, b):
    """
    merge forward and reverse drop-seq reads, annotating the reverse read (containing
    genomic information) with the cell barcode and rmt. Number of poly_t and pool fields
    are left empty.

    :param g: genomic fastq sequence data
    :param b: barcode fastq sequence data
    :return: annotated genomic sequence.
    """
    cell = b.sequence[:12]
    rmt = b.sequence[12:20]
    poly_t = b.sequence[20:-1]
    g.add_annotation((b'', cell, rmt, poly_t))
    return g


def mars1_seq(g, *args):
    """
    re-annotate reverse mars-seq reads in a format consistent with other SEQC platforms,
    annotating the reverse read (containing genomic information) with the pool, cell
    barcode, rmt, number of poly_t.

    :param g: genomic fastq sequence data
    :param b: barcode fastq sequence data
    :return: annotated genomic sequence.
    """

    *name_fields, pool, cell, rmt = g.name[1:-1].split(b':')
    g.name = (b'@' + b':'.join((pool, cell, rmt, b'')) + b';' +
              b':'.join(name_fields) + b'\n')
    return g


def mars2_seq(g, b):
    """
    re-annotate reverse mars-seq v2 reads in a format consistent with other SEQC
    platforms, annotating the reverse read (containing genomic information) with the pool,
    cell barcode, rmt, number of poly_t.

    :param g: genomic fastq sequence data
    :param b: barcode fastq sequence data
    :return: annotated genomic sequence.
    """
    pool = g.sequence.strip()[5:9]
    seq = b.sequence.strip()
    cell = seq[:7]
    rmt = seq[7:15]
    poly_t = seq[15:]
    g.add_annotation((b'', pool + cell, rmt, poly_t))
    return g
