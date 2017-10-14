import os
import shutil
import re
import numpy as np
import pandas as pd
import json
import phenograph
from matplotlib import pyplot as plt
from jinja2 import Environment, PackageLoader
from collections import OrderedDict, namedtuple
from weasyprint import HTML
from seqc.stats.tsne import TSNE
from sklearn.decomposition import PCA
from sklearn.linear_model import LinearRegression
from seqc import plot


ImageContent = namedtuple('ImageContent', ['image', 'caption', 'legend'])

TextContent = namedtuple('TextContent', ['text'])

DataContent = namedtuple('DataContent', ['keys', 'values'])


class Section:

    __slots__ = ['name', 'content', 'filename']

    def __init__(self, name, content, filename):
        """

        :param str name: Section name
        :param OrderedDict content: ordered dictionary containing keys corresponding to
          header information and values which are Content classes
        """
        self.name = name
        self.content = content
        self.filename = filename

    def render(self, prefix='', link_sections=None, index_section=None):
        """renders this sub-section

        :param str prefix: prefix for the rendered html file
        :param list link_sections: list of Section objects to link to in the sidebar
          of this section
        :param index_section: section which serves as index.html file
        :return None:
        """
        env = Environment(loader=PackageLoader('seqc.summary', 'templates'))
        section_template = env.get_template('section_content.html')
        if link_sections is None:
            link_sections = [self]
        rendered_section = section_template.render(
            sections=link_sections, section=self,
            index_section_link=index_section.filename)
        with open(prefix + self.filename, 'w') as f:
            f.write(rendered_section)

    @classmethod
    def from_alignment_summary(cls, alignment_summary, filename):
        """Create a summary section from an alignment summary file produced by STAR

        :param str alignment_summary: STAR alignment summary
        :param str filename: html file name for this section
        :return:
        """
        with open(alignment_summary, 'r') as f:
            data = f.read()
        categories = OrderedDict()

        def split_lines(block):
            keys, values = zip(*(l.split('|') for l in block.split('\n') if '|' in l))
            keys = [k.strip() for k in keys if k.strip()]
            values = [v.strip() for v in values if v.strip()]
            return DataContent(keys, values)

        # time and input read numbers
        time, _, data = data.partition('UNIQUE READS:')
        categories['Run Time'] = split_lines(time)

        # unique reads and splicing
        splicing, _, data = data.partition('MULTI-MAPPING READS:')
        categories['Unique Reads and Splicing'] = split_lines(splicing)

        # multimapping reads
        multimapping, _, unmapped = data.partition('UNMAPPED READS:')
        categories['Multimapping Reads'] = split_lines(multimapping)
        categories['Unmapped Reads'] = split_lines(unmapped)

        return cls('STAR Alignment Summary', categories, filename)

    @classmethod
    def from_status_filters(cls, ra, filename):
        """run after ReadArray is initialized and initial_filtering() has been run.

        :param ra: ReadArray object
        :param str filename: html file name for this section
        :return cls: Section containing initial filtering results
        """
        # todo replace whitespace characters with html equiv, add space b/w lines
        description = (
            'Initial filters are run over the sam file while our ReadArray database is '
            'being constructed. These filters indicate heuristic reasons why reads '
            'should be omitted from downstream operations:<br><br>'
            '<b>no gene</b>: Regardless of the read\'s genomic alignment status, there was no '
            'transcriptomic alignment for this read.<br>'
            '<b>gene not unique</b>: this indicates that more than one alignment was recovered '
            'for this read. We attempt to resolve these multi-alignments downstream. <br>'
            '<b>primer missing</b>: This is an in-drop specific filter, it indices that the '
            'spacer sequence could not be identified, and thus neither a cell barcode '
            'nor an rmt were recorded for this read.<br>'
            '<b>low poly t</b>: the primer did not display enough t-sequence in the primer '
            'tail, where these nucleotides are expected. This indicates an increased '
            'probability that this primer randomly primed, instead of hybridizing with '
            'the poly-a tail of an mRNA molecule.')
        description_section = TextContent(description)

        # Get counts 
        no_gene = np.sum(ra.data['status'] & ra.filter_codes['no_gene'] > 0)
        gene_not_unique = np.sum(ra.data['status'] & ra.filter_codes['gene_not_unique'] > 0)
        primer_missing = np.sum(ra.data['status'] & ra.filter_codes['primer_missing'] > 0)
        low_polyt = np.sum(ra.data['status'] & ra.filter_codes['low_polyt'] > 0)

        keys = ('length of read array',  'reads with genomic alignments alone', 'reads mapping to multiple genes', 'primer missing', 'low poly t')
        values = (
            len(ra.data),
            '%d (%.2f%%)' % (no_gene, no_gene / len(ra.data) * 100),
            '%d (%.2f%%)' % (gene_not_unique, gene_not_unique / len(ra.data) * 100),
            '%d (%.2f%%)' % (primer_missing, primer_missing / len(ra.data) * 100),
            '%d (%.2f%%)' % (low_polyt, low_polyt / len(ra.data) * 100),
        )
        data_section = DataContent(keys, values)
        return cls(
            'Initial Filtering',
            {'Description': description_section, 'Results': data_section},
            filename)

    @classmethod
    def from_cell_barcode_correction(cls, ra, filename):
        """Status page for cell barcode correction

        later, should add a figure for error rates, which will need to be returned by
        ra.apply_barcode_correction()

        :param ra:
        :param str filename: html file name for this section
        :return:
        """
        description = 'description for cell barcode correction'  # todo implement
        description_section = TextContent(description)
        count = np.sum(ra.data['status'] & ra.filter_codes['cell_error'] > 0)
        data_section = DataContent(
            ['cell error'],
            ['%d (%.2f%%)' % (count, count / len(ra.data) * 100)])
        return cls(
            'Cell Barcode Correction',
            {'Description': description_section, 'Results': data_section},
            filename)

    @classmethod
    def from_rmt_correction(cls, ra, filename):
        """Status page for error correction

        For now, returns the number of errors returned and a description of the rationale

        :param ra:
        :param str filename: html file name for this section
        :return:
        """

        description = 'description for rmt correction'  # todo implement
        description_section = TextContent(description)
        count = np.sum(ra.data['status'] & ra.filter_codes['rmt_error'])
        data_section = DataContent(
            ['rmt error'],
            ['%d (%.2f%%)' % (count, count / len(ra.data) * 100)])
        return cls(
            'RMT Barcode Correction',
            {'Description': description_section, 'Results': data_section},
            filename)

    @classmethod
    def from_resolve_multiple_alignments(cls, results, filename):
        """

        reports the number of corrected alignments

        :param dict results:
        :param str filename: html file name for this section
        :return:
        """
        description = 'description for multialignment correction'  # todo implement
        description_section = TextContent(description)
        keys, values = zip(*results.items())
        data_section = DataContent(keys, values)
        return cls(
            'Multialignment Resolution',
            {'Description': description_section, 'Results': data_section},
            filename)

    @classmethod
    def from_cell_filtering(cls, figure_path, filename):
        """

        This is the figure, but it needs a caption!

        :param str figure_path:
        :param str filename: html file name for this section
        :return:
        """
        description = 'description for cell filtering'  # todo implement
        description_section = TextContent(description)
        image_legend = 'image legend'  # todo implement
        image_section = ImageContent(figure_path, 'cell filtering figure', image_legend)
        return cls(
            'Cell Filtering',
            {'Description': description_section, 'Results': image_section},
            filename)

    @classmethod
    def from_final_matrix(cls, counts_matrix, figure_path, filename):
        """Create a histogram of cell sizes, a tSNE projection, and a diffusion map.

        save this data in a .h5 file accessible to pandas.

        :param filename:
        :param figure_path:
        :param pd.DataFrame counts_matrix:
        :return:
        """
        plot.Diagnostics.cell_size_histogram(counts_matrix, save=figure_path)
        
        # Number of cells and molecule count distributions
        image_legend = "Number of cells: {} <br>".format(counts_matrix.shape[0])
        ms = counts_matrix.sum(axis=1)
        image_legend += "Min number of molecules: {}<br>".format(ms.min())
        for prctile in [25, 50, 75]:
            image_legend += '{}th percentile: {}<br>'.format(prctile, np.percentile(ms, prctile))
        image_legend += "Max number of molecules: {}<br>".format(ms.max())

        image_section = ImageContent(figure_path, 'cell size figure', image_legend)
        return cls('Cell Summary',
                   {'Library Size Distribution': image_section},
                   filename)

    @classmethod
    def from_run_time(cls, log_file, filename):
        """

        :param filename:
        :param log_file:  seqc log file
        :return:
        """
        with open(log_file) as f:
            log = f.readlines()
        text_section = TextContent('<br>'.join(log))
        return cls('SEQC Log', {'Log Content': text_section}, filename)

    @classmethod
    def from_basic_clustering_and_projection(cls):
        """What if anything do we want to do here?

        # parameterize (default = True?)
        user could specify desire for median norm, PCA + tsne.

        :return:
        """
        raise NotImplementedError

    @classmethod
    def from_overall_yield(cls):
        """

        % losses from original yield summary, can wait.

        :return:
        """
        raise NotImplementedError


class Summary:

    def __init__(self, archive_name, sections, index_section=None):
        """
        :param str archive_name: filepath for the archive to be constructed
        :param list sections: dictionary of str filename: Section objects
        :param index_section: section to be produced as the index.html page.
        """
        self.archive_name = archive_name
        self.sections = sections
        self.reference_directory = os.path.dirname(__file__)
        self.index_section = index_section

    def prepare_archive(self):
        """

        :return:
        """
        if os.path.isdir(self.archive_name):
            shutil.rmtree(self.archive_name)
        shutil.copytree(self.reference_directory, self.archive_name)

    def import_image(self, image_path):
        filename = os.path.split(image_path)[1]
        shutil.copy(image_path, self.archive_name + '/img/' + filename)

    def render(self):
        """loop over sections and render them to filename keys within archive_name

        :return str zipped_archive_location:
        """
        html_location = self.archive_name + '/html_/'
        self.index_section.render(html_location, self.sections, self.index_section)
        for section in self.sections:
            section.render(html_location, self.sections, self.index_section)
        # todo this is not currently working.
        # os.symlink(html_location + self.index_section.filename,
        #            self.archive_name + '/index.html')

    def compress_archive(self):
        root_dir, _, base_dir = self.archive_name.rpartition('/')
        shutil.make_archive(
            self.archive_name, 'gztar', root_dir, base_dir)
        return self.archive_name + '.tar.gz'


class MiniSummary:
    def __init__(self, output_prefix, mini_summary_d, alignment_summary_file, filter_fig, cellsize_fig):
        """
        :param mini_summary_d: dictionary containing output parameters
        :param count_mat: count matrix after filtered
        :param filter_fig: filtering figure
        :param cellsize_fig: cell size figure
        """
        self.output_prefix = output_prefix
        self.mini_summary_d = mini_summary_d
        self.alignment_summary_file = alignment_summary_file
        self.filter_fig = filter_fig
        self.cellsize_fig = cellsize_fig
        self.pca_fig = output_prefix+"_pca.png"
        self.tsne_and_phenograph_fig = output_prefix+"_phenograph.png" 

    def compute_summary_fields(self, read_array, count_mat):
        self.count_mat = pd.DataFrame(count_mat)
        self.mini_summary_d['unmapped_pct'] = 0.0
        with open(self.alignment_summary_file, "r") as f:
            for line in f:
                arr = re.split("[\|:\t ]+", line)
                if "Number of input reads" in line:
                    self.mini_summary_d['n_reads'] = int(arr[-1].strip())
                elif "Uniquely mapped reads %" in line:
                    self.mini_summary_d['uniqmapped_pct'] = float(arr[-1].strip().strip("%"))
                elif "% of reads mapped to multiple loci" in line:
                    self.mini_summary_d['multimapped_pct'] = float(arr[-1].strip().strip("%"))
                elif "% of reads unmapped: too many mismatches" in line:
                    self.mini_summary_d['unmapped_pct'] += float(arr[-1].strip().strip("%"))
                elif "% of reads unmapped: too short" in line:
                    self.mini_summary_d['unmapped_pct'] += float(arr[-1].strip().strip("%"))
                elif "% of reads unmapped: other" in line:
                    self.mini_summary_d['unmapped_pct'] += float(arr[-1].strip().strip("%"))

        no_gene = np.sum(read_array.data['status'] & read_array.filter_codes['no_gene'] > 0)
        self.mini_summary_d['genomic_read_pct'] = no_gene / len(read_array.data) * 100

        # Calculate statistics from count matrix
        self.mini_summary_d['med_molcs_per_cell'] = np.median(count_mat.sum(1))
        self.mini_summary_d['molcs_per_cell_25p'] = np.percentile(count_mat.sum(1), 25)
        self.mini_summary_d['molcs_per_cell_75p'] = np.percentile(count_mat.sum(1), 75)
        self.mini_summary_d['molcs_per_cell_min'] = np.asscalar(np.min(count_mat.sum(1)))
        self.mini_summary_d['molcs_per_cell_max'] = np.asscalar(np.max(count_mat.sum(1)))
        self.mini_summary_d['n_cells'] = len(count_mat.index)

        # Filter low occurrence genes and median normalization
        self.counts_filtered = self.count_mat.loc[:, (self.count_mat>0).sum(0) >= 10]
        median_counts = np.median(self.counts_filtered.sum(1))
        counts_normalized = self.counts_filtered.divide(self.counts_filtered.sum(1),axis=0).multiply(median_counts)

        # Doing PCA transformation
        pcaModel = PCA(n_components=20)
        counts_pca_reduced = pcaModel.fit_transform(counts_normalized.as_matrix())

        # taking at most 20 components or total variance is greater than 80%
        num_comps = 0
        total_variance = 0.0
        while (total_variance < 80.0) and (num_comps < len(pcaModel.explained_variance_ratio_)):
            total_variance += pcaModel.explained_variance_ratio_[num_comps]
            num_comps += 1

        self.counts_after_pca = counts_pca_reduced[:, :num_comps]
        self.explained_variance_ratio = pcaModel.explained_variance_ratio_

        # regressed library size out of principal components 
        for c in range(num_comps):
            lm = LinearRegression(normalize=False)
            X = self.counts_filtered.sum(1).values.reshape(len(self.counts_filtered), 1)
            Y = counts_pca_reduced[:, c]
            lm.fit(X, Y)
            if c == 0:
                self.counts_pca_regressed_out_lib_size = Y-lm.predict(X)
            else:
                self.counts_pca_regressed_out_lib_size = np.column_stack((self.counts_pca_regressed_out_lib_size,
                                                                         Y - lm.predict(X)))

        # Doing TSNE transformation
        tsne = TSNE(n_components=2)
        self.counts_after_tsne = tsne.fit_transform(self.counts_pca_regressed_out_lib_size)
        self.clustering_communities, _, _ = phenograph.cluster(self.counts_pca_regressed_out_lib_size, k=50)

    def render(self):
        plot.Diagnostics.pca_components(self.pca_fig, self.explained_variance_ratio, self.counts_after_pca)
        plot.Diagnostics.phenograph_clustering(self.tsne_and_phenograph_fig, self.counts_filtered.sum(1), 
                                               self.clustering_communities, self.counts_after_tsne)

        self.mini_summary_d['seq_sat_rate'] = ((self.mini_summary_d['avg_reads_per_molc'] - 1.0) * 100.0
                                              / self.mini_summary_d['avg_reads_per_molc'])

        warning_d = dict()
        if ('mt_rna_fraction' in self.mini_summary_d) and (self.mini_summary_d['mt_rna_fraction'] >= 30):
            warning_d["High percentage of cell death"] = "Yes (%.2f%%)" \
                                                        % (self.mini_summary_d['mt_rna_fraction'])
        else:
            warning_d["High percentage of cell death"] = "No"
        warning_d["Noisy first few principle components"] = "Yes" if (self.explained_variance_ratio[0]<=0.05) else "No"
        if self.mini_summary_d['seq_sat_rate'] <= 5.00:
            warning_d["Low sequencing saturation rate"] = ("Yes (%.2f%%)" % (self.mini_summary_d['seq_sat_rate'])) 
        else:
            warning_d["Low sequencing saturation rate"] = "No"

        env = Environment(loader=PackageLoader('seqc.summary', 'templates'))
        section_template = env.get_template('mini_summary_base.html')
        rendered_section = section_template.render(output_prefix = self.output_prefix, warning_d = warning_d, 
                                                   mini_summary_d = self.mini_summary_d, cellsize_fig = self.cellsize_fig,
                                                   pca_fig = self.pca_fig, filter_fig = self.filter_fig,
                                                   tsne_and_phenograph_fig = self.tsne_and_phenograph_fig)
        with open(self.output_prefix + "_mini_summary.html", 'w') as f:
            f.write(rendered_section)

        HTML(self.output_prefix + "_mini_summary.html").write_pdf(self.output_prefix + "_mini_summary.pdf")

        with open(self.output_prefix + "_mini_summary.json","w") as f:
            json.dump(self.mini_summary_d, f)
        return self.output_prefix + "_mini_summary.json", self.output_prefix + "_mini_summary.pdf"