import numpy as np
import pandas as pd
from seqc.alignment import sam
from seqc.sequence.encodings import DNA3Bit
from scipy.sparse import coo_matrix, csr_matrix
import seqc.sequence.barcodes
import tables as tb
from itertools import permutations
from seqc import multialignment
from seqc.sparse_frame import SparseFrame
from seqc import log
from scipy.stats import hypergeom
from collections import OrderedDict


class ReadArray:

    _dtype = [
        ('status', np.uint8),  # if > 8 tests, change to int16
        ('cell', np.int64),
        ('rmt', np.int32),
        ('n_poly_t', np.uint8)]

    def __init__(self, data, genes, positions, chromosomes):
        """
        Enhanced np.ndarray (structured array) with several companion functions to hold
        filter, sort, and access compressed fastq read information.

        :param np.recarray data: np.ndarray with dtype of self._dtype
        :param csr_matrix|np.ndarray genes:
        :param csr_matrix|np.ndarray positions:

        """

        if not isinstance(genes, (csr_matrix, np.ndarray)):
            raise TypeError(
                'genes must be a scipy csr_matrix or np.array, not %s'
                % repr(type(genes)))
        self._genes = genes
        if not isinstance(positions, (csr_matrix, np.ndarray)):
            raise TypeError(
                'positions must be a scipy csr_matrix or np.array, not %s'
                % repr(type(positions)))
        self._positions = positions
        if not isinstance(chromosomes, (csr_matrix, np.ndarray)):
            raise TypeError(
                'chromosomes must be a scipy csr_matrix or np.array, not %s'
                % repr(type(chromosomes)))
        self._chromosomes = chromosomes
        if not isinstance(data, np.ndarray):
            raise TypeError(
                'data must be a structured np.array object, not %s' % repr(type(data)))
        self._data = data

        if isinstance(genes, csr_matrix):  # track if genes/positions are csr or np.array
            self._ambiguous_genes = True
        else:
            self._ambiguous_genes = False

    @property
    def data(self):
        return self._data

    @property
    def genes(self):
        return self._genes

    @genes.setter
    def genes(self, item):
        # todo: error check
        self._genes = item

    @property
    def positions(self):
        return self._positions

    @positions.setter
    def positions(self, item):
        # todo: error check
        self._positions = item
    
    @property
    def chromosomes(self):
        return self._chromosomes

    @chromosomes.setter
    def chromosomes(self, item):
        #todo: error check
        self._chromosomes = item

    def __len__(self):
        return len(self.data)

    def __iter__(self):
        """Returns an interator over the RecordArray

        If gene or position are zero, this indicates that those fields are empty. If the
        ReadArray has not yet been disambiguated, this function returns the first gene and
        position in each multialignment, if one is present.

        :return Iterator: iterator over (recarray_row, gene, position)
        """
        # todo test out speed loss with chain(), compare to numpy.concatenate, None
        if self._ambiguous_genes:  # genes, positions are scipy.sparse.csr_matrix
            for i in np.arange(self.data.shape[0]):
                yield self.data[i], self.genes[i, 0], self.positions[i, 0]
        else:  # genes, positions are np.array
            for i in np.arange(self.data.shape[0]):
                yield self.data[i], self.genes[i], self.positions[i]

    filter_codes = {
        'no_gene': 0b1,
        'rmt_error': 0b10,
        'cell_error': 0b100,  # todo this must execute before multialignment
        'low_polyt': 0b1000,
        'gene_not_unique': 0b10000,
        'primer_missing': 0b100000,
        'lonely_triplet': 0b1000000,  # todo could call this low coverage?
    }

    def initial_filtering(self, required_poly_t=1):
        """Apply different filters to the read array. If a read fails a filter, it is not
        passed to the others and so the counts for each filter is the number of reads that
        failed that one but passed all previous ones Filters are not ordered in any
        particular way.

        :param required_poly_t: minimum number of T nucleotides in the primer tail of a
          valid read
        :return None: method sets the status vector according to the result of filtering.
        """

        failing = np.zeros(self.data.shape[0], dtype=np.int8)

        # genes are dealt with differently depending on the state of the array
        if self._ambiguous_genes:
            nnz = self.genes.getnnz(axis=1)
            failing[nnz == 0] |= self.filter_codes['no_gene']
            failing[nnz > 1] |= self.filter_codes['gene_not_unique']
        else:  # multiple gene filter is empty
            failing[self.genes == 0] |= self.filter_codes['no_gene']

        # todo add logic for "primer_missing"
        failing[self.data['rmt'] == 0] |= self.filter_codes['primer_missing']
        failing[self.data['cell'] == 0] |= self.filter_codes['primer_missing']
        failing[self.data['n_poly_t'] < required_poly_t] |= self.filter_codes['low_polyt']

        self.data['status'] = np.bitwise_or(self.data['status'], failing)

    def filtering_mask(self, *ignore):
        """return a filtering mask that, when compared to status, will return False for
        reads failing any filters that are not specified in ignore

        :param *str ignore: fields to ignore when constructing filtering mask
        :return int: mask
        """
        mask = 2 ** 8 - 1  # 0b111... x16, permissive mask
        for filter_ in ignore:
            try:
                mask ^= self.filter_codes[filter_]  # mask filter
            except KeyError:
                raise KeyError('%s is not a valid filter. Please select from %s' %
                               (filter_, repr(self.filter_codes.keys())))
        return mask

    def iter_active(self, *ignore):
        """Iterator over active reads, ignoring any filter passed in ignore

        :param *str ignore: values to ignore when parsing active reads. choices:
          [no_gene, no_rmt, no_cell, low_polyt, gene_not_unique, no_spacer]
        :yields int, (np.array, int, int): iterator yields active records in form
          (index, (data_row, gene, pos))
        """

        if not ignore:  # save a bit of work by not &ing with mask
            for i, (data, gene, position, chromosome) in enumerate(self):
                if data['status'] == 0:
                    yield i, data, gene, position, chromosome

        else:  # create the appropriate mask for filters we want to ignore
            mask = self.filtering_mask(*ignore)
            for i, (data, gene, position, chromosome) in enumerate(self):
                if not data['status'] & mask:  # ignores fields in ignore
                    yield i, data, gene, position, chromosome

    @classmethod
    def from_alignment_file(cls, alignment_file, translator, required_poly_t):
        """
        construct a ReadArray object from a samfile containing only uniquely aligned
        records

        :param required_poly_t: number of poly_t required for a read to be considered
          a valid alignment
        :param GeneIntervals translator: translator created from the .gtf annotation
          file corresponding to the genome against which the reads in sam_file were
          aligned
        :param str alignment_file: filename of alignment file.
        :return:
        """

        # todo add a check for @GO query header (file matches sorting assumptions)

        reader = sam.Reader(alignment_file)  # todo swap to pysam reader, probably faster

        # todo allow reading of this from alignment summary
        num_reads = 0
        num_unique = 0
        prev_alignment_name = ''
        for alignment in reader:
            num_reads += 1
            if alignment.qname != prev_alignment_name:
                num_unique += 1
                prev_alignment_name = alignment.qname

        # pre-allocate arrays
        data = np.recarray((num_unique,), cls._dtype)
        row = np.zeros(num_reads, dtype=np.int32)
        col = np.zeros(num_reads, dtype=np.int8)
        position = np.zeros(num_reads, dtype=np.int32)
        gene = np.zeros(num_reads, dtype=np.int32)
        chromosome = np.zeros(num_reads, dtype=np.int32)

        # loop over multialignments
        row_idx = 0  # identifies the read index
        arr_idx = 0  # identifies the alignment index across all reads
        max_ma = 0  # maximum number of alignments observed for a read
        for ma in reader.iter_multialignments():
            col_idx = 0  # identifies the alignment index within a read
            for a in ma:
                genes = translator.translate(a.rname, a.strand, a.pos)

                # if gene passes filter, store passing genes and increment appropriate
                # arr and col indices
                if genes is not None:
                    row[arr_idx] = row_idx
                    col[arr_idx] = col_idx
                    gene[arr_idx] = genes
                    position[arr_idx] = a.pos
                    chromosome[arr_idx] = a.rname
                    arr_idx += 1
                    col_idx += 1
            max_ma = max(max_ma, col_idx)

            cell = seqc.sequence.encodings.DNA3Bit.encode(a.cell)
            rmt = seqc.sequence.encodings.DNA3Bit.encode(a.rmt)
            n_poly_t = a.poly_t.count('T') + a.poly_t.count('N')
            data[row_idx] = (0, cell, rmt, n_poly_t)
            row_idx += 1

        # some reads will not have aligned, throw away excess allocated space before
        # creating the ReadArray
        row, col, position, chromosome, gene = (
            row[:arr_idx], col[:arr_idx], position[:arr_idx], chromosome[:arr_idx], gene[:arr_idx])

        gene = coo_matrix((gene, (row, col)), shape=(row_idx, max_ma), dtype=np.int32)
        position = coo_matrix((position, (row, col)), shape=(row_idx, max_ma),
                              dtype=np.int32)
        chromosome = coo_matrix((chromosome, (row, col)), shape=(row_idx, max_ma),
                                dtype=np.int32)

        ra = cls(data, gene.tocsr(), position.tocsr(), chromosome.tocsr())
        ra.initial_filtering(required_poly_t=required_poly_t)
        return ra

    def group_indices_by_cell(self, multimapping=False):
        """group the reads in ra.data by cell.

        :param bool multimapping: if True, then reads are not filtered if they have
          alignments to more than one gene
        :return [np.array]: list of numpy arrays, each containing all of the indices for
          reads that correspond to a group, defined as a unique combination of the
          columns specified in parameter by.
        """
        idx = np.argsort(self.data['cell'])

        # filter the index for reads that
        if multimapping:
            mask = self.filtering_mask('gene_not_unique')
            passing = (self.data['status'][idx] & mask) == 0
        else:
            passing = self.data['status'][idx] == 0
        idx = idx[passing]

        # determine which positions in idx are the start of new groups (boolean, True)
        # convert boolean positions to indices, add start and end points.
        breaks = np.where(np.diff(self.data['cell'][idx]))[0] + 1
        # use these break points to split the filtered index according to "by"
        return np.split(idx, breaks)


    def save(self, archive_name):
        """save a ReadArray object as an hdf5 archive

        :param str archive_name: filestem for the new archive
        :return None:
        """

        def store_carray(archive, array, name):
            atom = tb.Atom.from_dtype(array.dtype)
            store = archive.create_carray(archive.root, name, atom, array.shape)
            store[:] = array
            store.flush()

        if not archive_name.endswith('.h5'):
            archive_name += '.h5'

        # construct container
        blosc5 = tb.Filters(complevel=5, complib='blosc')
        f = tb.open_file(archive_name, mode='w', title='Data for seqc.ReadArray',
                         filters=blosc5)

        f.create_table(f.root, 'data', self.data)

        if self._ambiguous_genes:
            # each array is data, indices, indptr
            store_carray(f, self.genes.indices, 'indices')
            store_carray(f, self.genes.indptr, 'indptr')
            store_carray(f, self.genes.data, 'gene_data')
            store_carray(f, self.positions.data, 'positions_data')
            store_carray(f, self.chromosomes.data, 'chromosome_data')
        else:
            store_carray(f, self.genes, 'genes')
            store_carray(f, self.positions, 'positions')
            store_carray(f, self.chromosomes, 'chromosomes')

        f.close()

    @classmethod
    def load(cls, archive_name):
        """load a ReadArray from an hdf5 archive, note that ma_pos and ma_genes are
        discarded.

        :param str archive_name: name of a .h5 archive containing a saved ReadArray object
        :return ReadArray:
        """

        f = tb.open_file(archive_name, mode='r')
        data = f.root.data.read()

        try:
            f.get_node('/genes')
            genes = f.root.genes.read()
            positions = f.root.positions.read()
            chromosomes = f.root.chromosomes.read()
        except tb.NoSuchNodeError:
            indptr = f.root.indptr.read()
            indices = f.root.indices.read()
            genes = f.root.gene_data.read()
            positions = f.root.positions_data.read()
            chromosomes = f.root.chromosome_data.read()
            genes = csr_matrix((genes, indices, indptr))
            positions = csr_matrix((positions, indices, indptr))
            chromosomes = csr_matrix((chromosomes, indices, indptr))

        return cls(data, genes, positions, chromosomes)

    # todo document me
    def resolve_ambiguous_alignments(self):
        """

        :return dict mm_results: dictionary containing information on how many molecules
          were resolved by this algorithm
        """

        # Resolve alignments
        indices_grouped_by_cells = self.group_indices_by_cell(multimapping=True)
        mm_results = self._resolve_alignments(indices_grouped_by_cells)

        # Reset genes and positions to be an array
        self._ambiguous_genes = False
        self.genes = np.ravel(self.genes.tocsc()[:, 0].todense()) 
        self.positions = np.ravel(self.positions.tocsc()[:, 0].todense())
        self.chromosomes = np.ravel(self.chromosomes.tocsc()[:, 0].todense())
        return mm_results

    # todo reconcile documentation with new method
    def _resolve_alignments(self, indices_grouped_by_cells):
        """
        Resolve ambiguously aligned molecules and edit the ReadArray data structures
        in-place to reflect the more specific gene assignments.

        After loading the co alignment matrix we group the reads of the ra by cell/rmt. 
        In each group we look at the different disjoint subsetes of genes reads are
        aligned to.

        side effect: # todo write here.

        :param list indices_grouped_by_cells: list of numpy arrays containing indices to
          ReadArray rows that correspond to each cell
        :return dict results: dictionary containing information on how many molecules
          were resolved by this algorithm
        """
        
        # Mask for reseting status on resolved genes
        mask = self.filtering_mask('gene_not_unique')

        # results dictionary for tracking effect of algorithm
        results = OrderedDict((
            ('unique molecules', 0),
            ('cell/rmt barcode collisions', 0),
            ('resolved molecules: disjoint', 0),
            ('resolved molecules: model', 0),
            ('ambiguous molecules', 0)
        ))

        for cell_group in indices_grouped_by_cells:

            # Sort by molecules
            inds = cell_group[np.argsort(self.data['rmt'][cell_group])]
            breaks = np.where(np.diff(self.data['rmt'][inds]))[0] + 1
            indices_grouped_by_molecule = np.split(inds, breaks)

            # Each molecule group
            for mol_group in indices_grouped_by_molecule:

                # Number of reads mapping to different sets of genes
                gene_groups = {}
                for ind in mol_group:
                    key = tuple(np.sort(self.genes[ind].data))
                    try:
                        gene_groups[key].append(ind)
                    except KeyError:
                        gene_groups[key] = [ind]

                # Return if there is only one gene group
                if len(gene_groups) == 1:
                    results['unique molecules'] += 1
                    continue

                # if it was not unique, there is a collision
                results['cell/rmt barcode collisions'] += 1

                # Divide into disjoint sets
                uf = multialignment.UnionFind()
                uf.union_all(gene_groups.keys())
                set_membership, sets = uf.find_all(gene_groups.keys())

                # Disambiguate each set
                keys = np.array(list(gene_groups.keys()))
                for s in sets:
                    set_groups = keys[set_membership == s]

                    # Return if the set contains only one group
                    if len(set_groups) == 1:
                        results['resolved molecules: disjoint'] += 1
                        continue

                    # Disambiguate if possible
                    common = list(multialignment.intersection(set_groups))

                    # Resolved if there is only one common gene
                    if len(common) == 1:
                        results['resolved molecules: model'] += 1
                        for group in set_groups:
                            # Update status
                            self.data['status'][gene_groups[tuple(group)]] &= mask
                            for ind in gene_groups[tuple(group)]:
                                # Update gene and position
                                if self.genes[ind, 0] == common[0]:
                                    continue
                                gene_index = (
                                    self.genes[ind] == common[0]).nonzero()[1][0]
                                self.genes[ind, 0] = self.genes[ind, gene_index]
                                self.positions[ind, 0] = self.positions[ind, gene_index]
                                self.chromosomes[ind, 0] = self.chromosomes[ind, gene_index]
                    else:
                        results['ambiguous molecules'] += 1
                        # Todo: Likelihood model goes here
        return results

    # todo : document me
    # Triplet filter from Adam
    def filter_low_coverage(self, alpha=0.25):
        
        use_inds = np.where( self.data['status'] == 0 )[0]
        cell = self.data['cell'][use_inds]
        position = self.positions[use_inds]
        chromosome = self.chromosomes[use_inds]
        rmt = self.data['rmt'][use_inds]
        genes = self.genes[use_inds]
        
        # A triplet is a (cell, position, rmt) triplet in each gene
        df = pd.DataFrame({'gene': genes, 'cell': cell, 'position': position, 
                           'rmt': rmt, 'chromosome': chromosome})
        grouped = df.groupby(['gene', 'position'])
        # This gives the gene followed by the number of triplets at each position
        # Summing across each gene will give the number of total triplets in gene
        num_per_position = (grouped['position'].agg({
            'Num Triplets at Pos': np.count_nonzero})).reset_index() 
        
        
        # Total triplets in each gene
        trips_in_gene = (num_per_position.groupby(['gene'])
            )['Num Triplets at Pos'].agg({'Num Triplets at Gene': np.sum})
            
        trips_in_gene = trips_in_gene.reset_index()
        
        num_per_position = num_per_position.merge(trips_in_gene,how = 'left')
        
        
        # for each (c,rmt) in df check in grouped2 if it is lonely 
        # determine number of lonely triplets at each position    
        grouped2 = df.groupby(['gene','cell','rmt'])    
        # lonely_triplets = grouped2["position"].apply(lambda x: len(x.unique()))     
        # This is a list of each gene, cell, rmt combo and the positions with that criteria
        lonely_triplets = grouped2['position'].apply(np.unique)
        lonely_triplets = pd.DataFrame(lonely_triplets)
        
        # if the length is one, this is a lonely triplet
        lonely_triplets_u = lonely_triplets['position'].apply(len) 
        lonely_triplets_u = pd.DataFrame(lonely_triplets_u)
        
        lonely_triplets_u = lonely_triplets_u.reset_index()
        lonely_triplets = lonely_triplets.reset_index()
            
        # Rename the columns
        lonely_triplets = lonely_triplets.rename(columns=lambda x: x.replace(
            'position', 'lonely position'))    
        lonely_triplets_u = lonely_triplets_u.rename(columns=lambda x: x.replace(
            'position', 'num'))
        
        # merge the column that is the length of the positions array
        # take the ones with length 1
        lonely_triplets = lonely_triplets.merge(lonely_triplets_u,how = 'left')
        lonely_triplets = lonely_triplets.loc[lonely_triplets.loc[:,'num'] == 1,:]
        
        # This is the gene, cell, rmt combo and the position that is lonely
        # We need to convert the array to a scalar 
        scalar = lonely_triplets["lonely position"].apply(np.asscalar)
        lonely_triplets["lonely position"] = scalar
        # Now if we group as such, we can determine how many (c, rmt) paris exist at each position
        # This would be the number of lonely pairs at a position
        grouped3 = lonely_triplets.groupby(["gene","lonely position"])
        l_num_at_position = (grouped3["cell"].agg(['count'])).reset_index()
        l_num_at_position = l_num_at_position.rename(columns=lambda x: x.replace(
            'count', 'lonely triplets at pos'))    
        l_num_at_position = l_num_at_position.rename(columns=lambda x: x.replace(
            'lonely position', 'position'))
        # lonely pairs in each gene
        l_num_at_gene = (lonely_triplets.groupby(["gene"]))['lonely position'].agg(
            ['count'])
        l_num_at_gene = l_num_at_gene.reset_index()
        l_num_at_gene = l_num_at_gene.rename(columns=lambda x: x.replace(
            'count', 'lonely triplets at gen'))    

        # aggregate    
        total = l_num_at_position.merge(l_num_at_gene,how='left')
        total = total.merge(num_per_position, how = 'left')
               
        # scipy hypergeom
        p = total.apply(self._hypergeom_wrapper, axis=1)
        p = 1-p
        
        from statsmodels.sandbox.stats.multicomp import multipletests as mt
        adj_p = mt(p,alpha = alpha, method='fdr_bh')
        
        keep = pd.DataFrame(adj_p[0])
        total['remove'] = keep
        
        remove = total[total['remove'] == True]
        
        final = df.merge(remove, how="left")
        final = final[final["remove"] == True]

        # Indicies to remove
        remove_inds = use_inds[final.index.values]
        
        self.data['status'][remove_inds] |= self.filter_codes['lonely_triplet']


    def _hypergeom_wrapper(self, x):
            
        from scipy.stats import hypergeom
        p = hypergeom.cdf(x['lonely triplets at pos'],x['Num Triplets at Gene'],
                          x['lonely triplets at gen'],x['Num Triplets at Pos'])
        return p


    def to_count_matrix(self, csv_path=None, sparse_frame=False, genes_to_symbols=False):
        """Convert the ReadArray into a count matrix.

        Since the matrix is sparse we represent it with 3 columns: row (cell),
         col (gene), value if not 0 (reads/molecules)

        :param str csv_path: file prefix for .csv output
        :param bool sparse_frame: if True, return sparse frame objects
        :param genes_to_symbols: if SparseFrame is True, integer gene IDs are converted
          to symbols in the returned SparseFrame objects
        :return str, str, dict, dict:
          filename (read counts),
          filename (mol counts),
          dict of (cell, rmt) -> read counts
          dict of (cell, rmt) -> molecule counts
        """
        reads_mat = {}
        mols_mat = {}
        for i, data, gene, pos, chrom in self.iter_active():

            try:
                reads_mat[data['cell'], gene] += 1
            except KeyError:
                reads_mat[data['cell'], gene] = 1

            try:
                rmt = data['rmt']
                if rmt not in mols_mat[data['cell'], gene]:
                    mols_mat[data['cell'], gene].append(rmt)
            except KeyError:
                mols_mat[data['cell'], gene] = [rmt]

        if sparse_frame:
            return (SparseFrame.from_dict(reads_mat, genes_to_symbols=genes_to_symbols),
                    SparseFrame.from_dict(
                        {k: len(v) for k, v in mols_mat.items()},
                        genes_to_symbols=genes_to_symbols))

        if csv_path is None:
            return reads_mat, mols_mat

        # todo convert gene integers to symbols before saving csv
        f = open(csv_path+'reads_count.csv', 'w')
        for data['cell'], gene in reads_mat:
            f.write('{},{},{}\n'.format(data['cell'], gene, reads_mat[data['cell'], gene]))

        f.close()


    def to_density_matrix(self):
        """Convert the ReadArray into a position density matrix
        
        :return dict:
            dict of (cell, chrom) -> read counts
            #TODO add counts by position windows in addition to whole chromosome counts
        """

        
        reads_mat = {}

        for i, data, gene, pos, chrom in self.iter_active():

            try:
                reads_mat[data['cell'], chrom] += 1
            except KeyError:
                reads_mat[data['cell'] chrom] = 1

        return reads_mat
