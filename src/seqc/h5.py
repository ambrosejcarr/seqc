__author__ = 'ambrose'

import tables as tb


def create_group(h5file, group_name, where='/', title=None, createparents=True):
    try:
        return h5file.get_node(where, group_name)  # exists, do nothing.
    except tb.NoSuchNodeError:
        return h5file.create_group(where=where, name=group_name, title=title,
                                   createparents=createparents)
