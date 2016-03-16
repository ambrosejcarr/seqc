import regex as re


_pattern = re.compile(b'(.{8,11}?)(GAGTGATTGCTTGTGACGCCTT){s<=2}(.{8})(.{6})(.*)')
_pattern_v2 = re.compile(b'(.{8,11}?)(GAGTGATTGCTTGTGACGCCAA){s<=2}(.{8})(.{8})(.*)')


def _check_spacer_v2(b):
    """a fast, in-drop-v1 specific command to find a spacer sequence and cb length

    :b: barcode fastq
    :returns: (cell_barcode, rmt, poly_t)
    """
    identifier = b.sequence[28:32]
    if identifier == b'CGCC':
        cb1 = b.sequence[:8]
        cb2 = b.sequence[30:38]
        rmt = b.sequence[38:46]
        poly_t = b.sequence[46:]
    elif identifier == b'ACGC':
        cb1 = b.sequence[:9]
        cb2 = b.sequence[31:39]
        rmt = b.sequence[39:47]
        poly_t = b.sequence[47:]
    elif identifier == b'GACG':
        cb1 = b.sequence[:10]
        cb2 = b.sequence[32:40]
        rmt = b.sequence[40:48]
        poly_t = b.sequence[48:]
    elif identifier == b'TGAC':
        cb1 = b.sequence[:11]
        cb2 = b.sequence[33:41]
        rmt = b.sequence[41:49]
        poly_t = b.sequence[49:]
    else:
        return b'', b'', b''
    cell = cb1 + cb2
    return cell, rmt, poly_t


def _check_spacer_v1(b):
    """a fast, in-drop-v1 specific command to find a spacer sequence and cb length

    :arg b: barcode fastq
    :returns: (cell_barcode, rmt, poly_t)
    """
    identifier = b.sequence[28:32]
    if identifier == b'CGCC':
        cb1 = b.sequence[:8]
        cb2 = b.sequence[30:38]
        rmt = b.sequence[38:44]
        poly_t = b.sequence[44:]
    elif identifier == b'ACGC':
        cb1 = b.sequence[:9]
        cb2 = b.sequence[31:39]
        rmt = b.sequence[39:45]
        poly_t = b.sequence[45:]
    elif identifier == b'GACG':
        cb1 = b.sequence[:10]
        cb2 = b.sequence[32:40]
        rmt = b.sequence[40:46]
        poly_t = b.sequence[46:]
    elif identifier == b'TGAC':
        cb1 = b.sequence[:11]
        cb2 = b.sequence[33:41]
        rmt = b.sequence[41:51]
        poly_t = b.sequence[51:]
    else:
        return b'', b'', b''
    cell = cb1 + cb2
    return cell, rmt, poly_t


def in_drop(g, b):
    cell, rmt, poly_t = _check_spacer_v1(b)
    if not cell:
        try:
            cell1, spacer, cell2, rmt, poly_t = re.match(
                _pattern, b.sequence).groups()
            cell = cell1 + cell2
        except AttributeError:
            cell, rmt, poly_t = b'', b'', b''
    g.add_annotation((b'', cell, rmt, poly_t))
    return g


def in_drop_v2(g, b):
    """
    :param g:
    :param b:
    :return:
    """
    cell, rmt, poly_t = _check_spacer_v2(b)
    if not cell:
        try:
            cell1, spacer, cell2, rmt, poly_t = re.match(
                _pattern_v2, b.sequence).groups()
            cell = cell1 + cell2
        except AttributeError:
            cell, rmt, poly_t = b'', b'', b''
    g.add_annotation((b'', cell, rmt, poly_t))
    return g


def drop_seq(g, b):
    """
    :param g:
    :param b:
    :return:
    """
    cell = b.sequence[:12]
    rmt = b.sequence[12:20]
    poly_t = b.sequence[20:]
    g.add_annotation((b'', cell, rmt, poly_t))
    return g


def mars1_seq(g, *args):
    """
    :param g:
    :param args:
    :return:
    """
    *name_fields, pool, cell, rmt = g.name[1:-1].split(b':')
    g.name = (b'@' + b':'.join((pool, cell, rmt, b'')) + b';' +
              b':'.join(name_fields) + b'\n')
    return g


def mars2_seq(g, *args):
    """
    :param g:
    :param args:
    :return:
    """
    *name_fields, pool, cell, rmt = g.name[1:-1].split(b'-')
    g.name = (b'@' + b':'.join((pool, cell, rmt, b'')) + b';' +
              b':'.join(name_fields) + b'\n')
    return g
