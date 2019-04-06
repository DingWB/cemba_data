"""
Group visualization

Input:
- Tidy Data: row is cell/group, groupby columns (x, y), groupby orders (x, y)
- Matrix Data: row is x|y, col is y|x. groupby and groupby orders provide separately
- Dendrogram: provided separately

Plot Type:
- Dotplot/Matrixplot: tidy data, single axes, multi visual map (size, color, marker)
- Heatmap: matrix data, single ax, single visual map
- Stack-violin: matrix data, multi axes, single visual map
- TracePlot: matrix data, multi axes, single visual map

Appendix:
- dendrogram
- category color blocks: ax, orientation, order
- category brackets
"""
from .tree import dendrogram
import seaborn as sns


def fancy_heatmap(ax, tidy_data, row_col, col_col,
                  size_col, sig_cutoff, color_col, palette='viridis',
                  row_order=None, col_order=None,
                  heatmap_scatter_kws=None, sig_scatter_kws=None):
    data = tidy_data.copy()
    data[row_col] = data[row_col].astype('category')
    data[col_col] = data[col_col].astype('category')

    if row_order is None:
        row_dendro = dendrogram(data, groupby=row_col)
        row_order = row_dendro['dendrogram_info']['ivl']
    data['row_i'] = data[row_col].map(lambda i: row_order.index(i))

    if col_order is None:
        col_dendro = dendrogram(data, groupby=col_col, cor_method='spearman')
        col_order = col_dendro['dendrogram_info']['ivl']
    data['col_i'] = data[col_col].map(lambda i: col_order.index(i))
    data['sig_marker'] = data[size_col] > sig_cutoff

    # plot color blocks
    _heatmap_scatter_kws = dict(hue_norm=(-1, 1), size_norm=(sig_cutoff / 2, sig_cutoff * 2), sizes=(20, 800),
                                legend=None, marker='s')
    if heatmap_scatter_kws is not None:
        _heatmap_scatter_kws.update(heatmap_scatter_kws)
    sns.scatterplot(x='row_i', y='col_i',
                    data=data, palette=palette,
                    hue=color_col, size=size_col, **_heatmap_scatter_kws, ax=ax)
    # plot sig marker
    _sig_scatter_kws = dict(color='white', linewidth=1, s=100, marker='+')
    if sig_scatter_kws is not None:
        _sig_scatter_kws.update(sig_scatter_kws)
    sns.scatterplot(x='row_i', y='col_i',
                    data=data[data['sig_marker']], **_sig_scatter_kws, ax=ax)

    ax.set(xlim=(-0.5, 6.5),
           yticks=range(data[col_col].unique().size),
           yticklabels=col_order,
           xticks=range(data[row_col].unique().size),
           xticklabels=row_order)
    for label in ax.xaxis.get_ticklabels():
        label.set(rotation=45, rotation_mode='anchor', ha='right')
    sns.despine(ax=ax, left=True, bottom=True)
    return ax
