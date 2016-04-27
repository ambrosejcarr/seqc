import os
import warnings
import numpy as np
import matplotlib
from scipy.stats import gaussian_kde
try:
    os.environ['DISPLAY']
except KeyError:
    matplotlib.use('Agg')
import matplotlib.pyplot as plt
with warnings.catch_warnings():
    warnings.simplefilter('ignore')  # catch experimental ipython widget warning
    import seaborn as sns

# set plotting defaults
sns.set_style('ticks')
matplotlib.rc('font', **{
    'serif': ['Computer Modern Roman'],
    'monospace': ['Inconsolata'],
    'sans-serif': ['Lato']
    })

# todo change to kwargs
matplotlib.rc('figure', **{'figsize': (3, 3), 'dpi': 150})

matplotlib.rc('patch', **{'facecolor': 'royalblue', 'edgecolor': 'none'})

matplotlib.rc('lines', **{'color': 'royalblue', 'markersize': 7})

matplotlib.rc('savefig', **{'dpi': '150'})

matplotlib.rc('image', **{'cmap': 'viridis'})


class FigureGrid:

    def __init__(self, n: int, max_cols=3):
        self.n = n
        self.nrows = int(np.ceil(n / max_cols))
        self.ncols = int(min((max_cols, n)))
        figsize = self.ncols * 3, self.nrows * 3

        # create figure
        self.gs = plt.GridSpec(nrows=self.nrows, ncols=self.ncols)
        self.figure = plt.figure(figsize=figsize)

        # create axes
        self.axes = {}
        for i in range(n):
            row = int(i // self.ncols)
            col = int(i % self.ncols)
            self.axes[i] = plt.subplot(self.gs[row, col])

    def __getitem__(self, item):
        return self.axes[item]

    def __iter__(self):
        for i in range(self.n):
            yield self[i]

    def tight_layout(self):
        self.gs.tight_layout(self.figure)

    def despine(self, top=True, right=True, **kwargs):
        for i in range(self.n):
            sns.despine(ax=self[i], top=top, right=right, **kwargs)

    def savefig(self, *args, **kwargs):
        self.figure.savefig(*args, **kwargs)


def map_categorical_to_cmap(data: np.ndarray, cmap=plt.get_cmap()):
    """
    :param data:
    :param cmap:
    :return:
    """
    categories = np.unique(data)
    n = len(categories)
    colors = cmap(np.linspace(0, 1, n))
    category_to_color = dict(zip(categories, colors))
    return np.array([category_to_color[i] for i in data]), category_to_color


def add_legend_to_categorical_vector(
        colors: np.ndarray, labels, ax, loc='center left', bbox_to_anchor=(0.98, 0.5),
        markerscale=0.75, **kwargs):
    """
    :param colors:
    :param labels:
    :param ax:
    :param kwargs: additional kwargs for legend
    :return:
    """
    artists = []
    for c in colors:
        artists.append(plt.Line2D((0, 1), (0, 0), color=c, marker='o', linestyle=''))
    ax.legend(
        artists, labels, loc=loc, bbox_to_anchor=bbox_to_anchor, markerscale=markerscale,
        **kwargs)


class scatter:

    @staticmethod
    def categorical(
            x, y, c, ax=None, cmap=plt.get_cmap(), legend_kwargs=None, randomize=True,
            *args, **kwargs):
        """
        :param x:
        :param y:
        :param c:
        :param ax:
        :param cmap:
        :param legend_kwargs:
        :param args:
        :param kwargs:
        :return: ax
        """
        if not ax:  # todo replace with plt.gridspec() method
            ax = plt.gca()

        if legend_kwargs is None:
            legend_kwargs = dict()
        color_vector, category_to_color = map_categorical_to_cmap(c, cmap)

        if randomize:
            ind = np.random.permutation(len(x))
        else:
            ind = np.argsort(np.ravel(c))

        ax.scatter(np.ravel(x)[ind], np.ravel(y)[ind], c=color_vector[ind], *args,
                   **kwargs)

        labels, colors = zip(*sorted(category_to_color.items()))
        sns.despine(ax=ax)
        add_legend_to_categorical_vector(colors, labels, ax, **legend_kwargs)
        return ax

    @staticmethod
    def continuous(x, y, c=None, ax=None, colorbar=True, randomize=True, **kwargs):
        """
        :param x:
        :param y:
        :param c:
        :param ax:
        :param kwargs:
        :return: ax
        """

        if ax is None:
            ax = plt.gca()

        if c is None:  # plot density if no color vector is provided
            x, y, c = scatter.density_2d(x, y)

        if randomize:
            ind = np.random.permutation(len(x))
        else:
            ind = np.argsort(c)

        sm = ax.scatter(x[ind], y[ind], c=c[ind], **kwargs)
        ax.xaxis.set_major_locator(plt.NullLocator())
        ax.yaxis.set_major_locator(plt.NullLocator())
        if colorbar:
            cb = plt.colorbar(sm)
            cb.ax.xaxis.set_major_locator(plt.NullLocator())
            cb.ax.yaxis.set_major_locator(plt.NullLocator())
        return ax

    @staticmethod
    def density_2d(x, y):
        """return x and y and their density z, sorted by their density (smallest to largest)

        :param x:
        :param y:
        :return:
        """
        xy = np.vstack([np.ravel(x), np.ravel(y)])
        z = gaussian_kde(xy)(xy)
        return np.ravel(x), np.ravel(y), np.arcsinh(z)
