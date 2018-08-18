import h5py
import numpy as np
import pandas as pd
import json
import datetime
from sklearn.utils import sparsefuncs
from sklearn.preprocessing import Imputer
from ast import literal_eval
from scipy.sparse import csr_matrix, csc_matrix, issparse, lil_matrix, vstack, hstack
from anndata import AnnData, read_h5ad


def cur_time():
    return datetime.datetime.now().strftime("%Y%m%d-%H%M%S")


class Dataset:
    # TODO change the general sparse format into csc/csr. lil is no good
    # TODO learn the implementation of AnnData at this part first
    # TODO more about HDF5 links and the backed mode

    def __init__(self, h5mc_path, mode='r'):
        self.h5f = h5py.File(h5mc_path, mode=mode)
        return

    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc_value, traceback):
        self.h5f.__exit__(exc_type, exc_value, traceback)

    def close(self):
        self.h5f.close()

    def __getitem__(self, item):
        return self.h5f[item]

    @property
    def cell_meta(self):
        meta_df = pd.DataFrame(self.h5f['cell_meta'].value)
        meta_df.set_index('_id', inplace=True)
        return meta_df

    @property
    def assembly_args(self):
        arg_dict = json.loads(self.h5f.attrs['assembly_args'])
        eval_args = {}
        for k, v in arg_dict.items():
            try:
                eval_args[k] = literal_eval(v)
            except (SyntaxError, ValueError):
                eval_args[k] = v
        return eval_args

    def get_region_meta(self, region_name):
        meta_df = pd.DataFrame(self.h5f[region_name]['region_meta'].value)
        meta_df.set_index('_id', inplace=True)
        return meta_df

    # iter_cells_data and get_cells_matrix can be faster...
    def iter_cells_data(self, region_name, context, cells, sparse_format='coo'):
        context = context.upper()
        for cell in cells:
            query_path = f'/{region_name}/{context}/{cell}'
            cell_group = self.h5f[query_path]
            cell_mc = _parse_lil_group(cell_group['mc'], sparse_format)
            cell_cov = _parse_lil_group(cell_group['cov'], sparse_format)
            yield cell, cell_mc, cell_cov

    def get_cells_matrix(self, region_name, context, cells=None, sparse_format='csr', to_df=False):
        region_meta_df = self.get_region_meta(region_name)
        mc_list = []
        cov_list = []
        if cells is None:
            cells = self.cell_meta.index
        for cell, mc, cov in self.iter_cells_data(region_name, context, cells, sparse_format):
            mc_list.append(mc)
            cov_list.append(cov)
        mc_sparse_matrix = vstack(mc_list)  # concatenate by row
        cov_sparse_matrix = vstack(cov_list)
        if to_df:
            print(f'Return pandas DataFrame for {len(cells)} cells')
            mc_df = pd.DataFrame(mc_sparse_matrix.toarray(), columns=region_meta_df.index, index=cells)
            cov_df = pd.DataFrame(cov_sparse_matrix.toarray(), columns=region_meta_df.index, index=cells)
            return mc_df, cov_df
        else:
            column_index = region_meta_df.index
            row_index = cells
            return mc_sparse_matrix, cov_sparse_matrix, column_index, row_index

    # must used func
    def get_mc_rate(self, region_name, context, cov_cutoff, cells=None, study_name=None):
        mc_sparse_matrix, cov_sparse_matrix, column_index, row_index = \
            self.get_cells_matrix(region_name, context, cells)

        mask = cov_sparse_matrix > cov_cutoff
        cov_pass = cov_sparse_matrix.multiply(mask)  # value not pass cutoff will be 0
        mc_pass = mc_sparse_matrix.multiply(mask)  # use same filter as cov
        mc_rate = (mc_pass / cov_pass).astype(np.float32)
        # 0 is different from NA
        mc_rate[mc_rate == 0] = 1e-9  # assign 0 to a small number to be compatible with sparse matrix
        mc_rate[np.isnan(mc_rate)] = 0  # assign 0 to a small number to be compatible with sparse matrix

        row_dict = {'row_names': row_index}
        col_dict = {'col_names': column_index,
                    'context': [context] * len(column_index),
                    'cov_cutoff': [cov_cutoff] * len(column_index),
                    'region_set_name': [region_name] * len(column_index)}

        return Study(mc_rate=csr_matrix(mc_rate),
                     col_dict=col_dict,
                     row_dict=row_dict,
                     uns_dict=None,
                     study_name=study_name)


_COL_DICT_ESSENTIAL_KEYS = {'col_names',
                            'cov_cutoff',
                            'context',
                            'region_set_name'}
_ROW_DICT_ESSENTIAL_KEYS = {'row_names'}


class Study:
    """
    Study contains 4 major part:
    1. _mc_rate: mC%
    2. _row_dict: id and other information of rows
    3. _col_dict: id and other information of cols
    4. _uns_dict: general information

    For general analysis, currently I actually use AnnData and scanpy,
    because I don't think its necessary to rebuild the wheels.
    Study can be easily transferred into AnnData by .to_ann()

    The reason I write Study is to host some methods that are specific to methylation data.
    The output file of Study is actually AnnData too,
    which can also be load as a Study using prepare_study.read_from_ann()

    TODO Change col, row dict into dataframe
    TODO Distinguish inplace and copy clearly
    TODO Cooperate with the backed mode, currently everything should be in memory
    """

    def __init__(self, mc_rate=None, col_dict=None, row_dict=None, uns_dict=None, study_name=None):
        """
        If col_dict, row_dict, uns_dict provided, use these three directly, ignore corresponding args
        """
        # init from matrix and dicts
        # main data
        self._mc_rate = mc_rate
        # prepare col_dict
        # item check of col dict
        for k in _COL_DICT_ESSENTIAL_KEYS:
            if k not in col_dict:
                raise KeyError(f'{k} not found in col_dict.')
            else:
                if len(col_dict[k]) != mc_rate.shape[1]:
                    raise ValueError(f'Length of {k} in col_dict is not equal to the cols in data.')
        self._col_dict = col_dict
        # prepare row dict
        # item check of row dict
        for k in _ROW_DICT_ESSENTIAL_KEYS:
            if k not in row_dict:
                raise KeyError(f'{k} not found in row_dict.')
            else:
                if len(row_dict[k]) != mc_rate.shape[0]:
                    raise ValueError(f'Length of {k} in row_dict is not equal to the rows in data.')
        self._row_dict = row_dict
        if uns_dict is None:
            self._uns_dict = {'create_time': cur_time()}
            if study_name is not None:
                self._uns_dict['study_name'] = study_name
        else:
            self._uns_dict = uns_dict

        # for convenience
        self._row_idx = self._row_dict['row_names']
        self._col_idx = self._col_dict['col_names']

    def __add__(self, obj):
        """
        add should only be used in the first step, only take care of the essential attr,
        computed attr will lose (should because data changed)

        row (cell) concatenate of study.
        Three things need to be checked before concatenate rows:
        1. Region names
        2. Region mc_context
        3. Region cov_cutoff
        It make no sense adding cells if these three element is different.
        :param obj:
        :return:
        """
        if not isinstance(obj, Study):
            raise TypeError(f'Adding Study objects with {type(obj)} is not allowed.')
        # 1. Region names
        for x, y in zip(self._col_dict['col_names'], obj._col_dict['col_names']):
            if x != y:
                raise ValueError('Study objects must have same regions in same order')
        # 2. Region mc_context
        if set(self._col_dict['context']) != set(obj._col_dict['context']):
            raise ValueError('Adding Study objects with different context.')
        # 3. Region cov_cutoff
        if set(self._col_dict['cov_cutoff']) != set(obj._col_dict['cov_cutoff']):
            raise ValueError('Adding Study objects with different cov_cutoff')
        # col dict will not changed

        # row concatenate main data
        new_mc_rate = vstack([self._mc_rate.tocsr(), obj._mc_rate.tocsr()])
        # concatenate essential elements of row_dict
        # TODO change this into pd.concat
        new_row_dict = {k: np.concatenate([self._row_dict[k], obj._row_dict[k]])
                        for k in _ROW_DICT_ESSENTIAL_KEYS}

        # check new row idx, raise if duplicates found
        if len(set(new_row_dict['row_names'])) != len(new_row_dict['row_names']):
            raise ValueError('Concatenated cells have duplicate index.')

        # TODO add warning if row_dict or col_dict or uns_dict contain computed attr

        return Study(mc_rate=new_mc_rate, col_dict=self._col_dict, row_dict=new_row_dict, uns_dict=None)

    def __radd__(self, other):
        return self + other

    def __repr__(self):
        content = f'Study of {self.shape[0]} cells × {self.shape[1]} regions.\n' \
                  f'    Col features: {"; ".join(list(self._col_dict.keys()))}\n' \
                  f'        Region set include: {"; ".join(list(set(self._col_dict["region_set_name"])))}\n' \
                  f'        mC context include: {"; ".join(list(set(self._col_dict["context"])))}\n' \
                  f'        coverage cutoff include: {"; ".join(map(str, list(set(self._col_dict["cov_cutoff"]))))}\n' \
                  f'    Row features: {"; ".join(list(self._row_dict.keys()))}\n'
        return content

    def __getitem__(self, item):
        cls = type(self)

        row_dict = self._row_dict
        col_dict = self._col_dict
        uns_dict = self._uns_dict

        # the row and col dict select
        if isinstance(item, (slice, list)):
            row_dict = _slice_dict(self._row_dict, item)
            data = self._mc_rate[item, :]
        elif isinstance(item, tuple):
            # multiple axis
            if len(item) > 2:
                raise NotImplementedError('Index have {len(item)} dimensions, only accept 1 or 2')
            elif len(item) == 2:
                row_dict = _slice_dict(self._row_dict, item[0])
                col_dict = _slice_dict(self._col_dict, item[1])
                data = self._mc_rate.tocsr()[item[0], :]
                data = data.tocsc()[:, item[1]]
            else:
                return self.__getitem__(list(item))
        else:
            raise NotImplementedError('Index can not be int yet, should be a slice or list')
        # the data selection is direct passed to np/scipy
        return cls(mc_rate=data, col_dict=col_dict, row_dict=row_dict, uns_dict=uns_dict, study_name=None)

    def region_append(self, obj):
        """
        region_append should only be used in the first step, only take care of the essential attr,
        computed attr will lose (should because data changed)

        :param obj: MethylRate object
        :return:
        """
        # 1. check row
        if not isinstance(obj, Study):
            raise TypeError(f'region_append Study objects with {type(obj)}')
        for x, y in zip(self._row_idx, obj._row_idx):
            # strong check of index
            if x != y:
                raise ValueError('Study objects must have same cells in same order')

        # 2. concatenate col
        new_mc_rate = hstack([self._mc_rate.tocsc(), obj._mc_rate.tocsc()])
        # concatenate essential elements of col_dict
        new_col_dict = {k: np.concatenate([self._col_dict[k], obj._col_dict[k]])
                        for k in _COL_DICT_ESSENTIAL_KEYS}

        # check new col idx, warn if duplicates found
        if len(set(new_col_dict['col_names'])) != len(new_col_dict['col_names']):
            print('Warning: col_names have duplicates after region_append, '
                  'use make_col_name_unique() to modify that.')

        return Study(mc_rate=new_mc_rate, col_dict=new_col_dict, row_dict=self._row_dict, uns_dict=None)

    def to_ann(self):
        return AnnData(self._mc_rate, self._row_dict, self._col_dict, uns=self._uns_dict)

    def make_col_name_unique(self, use_key=None):
        if use_key is not None:
            # simply add number after dup col_name
            self._col_dict['col_names'] = self._col_dict['col_names']\
                                          + '_' + self._col_dict[use_key]
            self._col_idx = self._col_dict['col_names']
        if use_key is None:
            # add key items after dup col_name, if still have dup, warn and add number
            # TODO
            pass
        return

    def save(self, path, compression, compression_opts):
        ann = self.to_ann()
        ann.write(filename=path,
                  compression=compression,
                  compression_opts=compression_opts)
        return

    @classmethod
    def from_file(cls, path):
        ann = read_h5ad(path)  # load everything into memory
        col_dict = {i: c for i, c in ann.var.iteritems()}
        col_dict['col_names'] = ann.var.index
        row_dict = {i: c for i, c in ann.obs.iteritems()}
        row_dict['row_names'] = ann.obs.index
        return Study(mc_rate=ann.X, col_dict=col_dict, row_dict=row_dict, uns_dict=ann.uns)

    @property
    def value(self):
        return self._mc_rate

    @property
    def shape(self):
        return self._mc_rate.shape

    @property
    def cells(self):
        return self._row_idx

    @property
    def regions(self):
        return self._col_idx

    @property
    def rows(self):
        return self._row_dict

    @property
    def columns(self):
        return self._col_dict

    # Preprocess functions

    def filter_cell(self, max_na_rate):
        """
        add row_mask to row_dict under row_na_cutoff,
        mask is different from view,
        cell been masked mean to be ignored on any downstream analysis due to low quality.
        :param max_na_rate:
        :return:
        """
        self._uns_dict['row_na_cutoff'] = max_na_rate
        if 'col_mask' not in self._col_dict:
            self._row_dict['row_mask'], _stat = _get_sparse_na_mask(self._mc_rate, axis=1, cutoff=max_na_rate)
        else:
            col_mask = self._col_dict['col_mask']
            self._row_dict['row_mask'], _stat = _get_sparse_na_mask(self._mc_rate[:, np.ravel(col_mask)],
                                                                    axis=1, cutoff=max_na_rate)
        na_rate = '%.1f' % (sum(self._row_dict['row_mask']) / len(self._row_idx) * 100)
        print(f"{sum(self._row_dict['row_mask'])} cells "
              f"({na_rate}%, mean={_stat[0]:.2f}, std={_stat[1]:.2f}) remained after this cell filter")
        return

    def filter_region(self, max_na_rate):
        """
        add col_mask to col_dict under col_na_cutoff,
        mask is different from view,
        cell been masked mean to be ignored on any downstream analysis due to low quality.
        :param max_na_rate:
        :return:
        """
        self._uns_dict['col_na_cutoff'] = max_na_rate
        if 'row_mask' not in self._row_dict:
            self._col_dict['col_mask'], _stat = _get_sparse_na_mask(self._mc_rate, axis=0, cutoff=max_na_rate)
        else:
            row_mask = self._row_dict['row_mask']
            self._col_dict['col_mask'], _stat = _get_sparse_na_mask(self._mc_rate[np.ravel(row_mask), :],
                                                                    axis=0, cutoff=max_na_rate)
        na_rate = '%.1f' % (sum(self._col_dict['col_mask']) / len(self._col_idx) * 100)
        print(f"{sum(self._col_dict['col_mask'])} regions "
              f"({na_rate}%, mean={_stat[0]:.2f}, std={_stat[1]:.2f}) remained after this region filter")
        return

    def reset_filter(self):
        # use dict remove
        self._uns_dict.pop("col_na_cutoff", None)
        self._uns_dict.pop("row_na_cutoff", None)
        self._row_dict.pop("row_mask", None)
        self._col_dict.pop("col_mask", None)
        return

    def add_row_feature(self, name, data, allow_na=True):
        row_data = data.reindex(self._row_idx)
        if not allow_na:
            if row_data.isnull().sum() != 0:
                raise ValueError('NA in row feature.')
        self._row_dict[name] = data
        return

    def add_col_feature(self, name, data, allow_na=True):
        col_data = data.reindex(self._col_idx)
        if not allow_na:
            if col_data.isnull().sum() != 0:
                raise ValueError('NA in col feature.')
        self._col_dict[name] = data
        return

    def normalize_by_row_feature(self, feature_name):
        if feature_name not in self._row_dict:
            raise KeyError(f'{feature_name} not in row_dict')
        sparsefuncs.inplace_row_scale(self._mc_rate, 1 / self._row_dict[feature_name])
        print(f'mC rate is normalized on each row by {feature_name}, change happened inplace')
        return

    def scale(self):
        # TODO
        return

    def imputation(self, missing_values=0, strategy='mean', axis=0, verbose=0, copy=False):
        # Simple naive imputer
        Imputer(missing_values=missing_values,  # in mc_rate sparse matrix 0 is NaN, 1e-9 is 0
                strategy=strategy,
                axis=axis,  # 0, along column; 1, along row
                verbose=verbose,
                copy=copy  # inplace imputation
                ).fit_transform(self._mc_rate)
        print(f'Used sklearn.Imputer to impute by {strategy} along {"column" if axis == 0 else "row"}, '
              f'change happened inplace')
        return


def _get_sparse_na_mask(sparse_arr, axis, cutoff):
    """
    In the sparse matrix, real zero is set to 1e-9, na is set to 0, count 0 is count na
    :param sparse_arr:
    :param axis:
    :param cutoff:
    :return:
    """
    na_rate = 1 - (sparse_arr != 0).sum(axis=axis).A1 / sparse_arr.shape[axis]
    na_stat = (np.mean(na_rate), np.std(na_rate))
    na_mask = na_rate < cutoff
    return na_mask, na_stat


def _parse_lil_group(lil_group, sparse_format):
    data = lil_group['lil_data'].value
    index = lil_group['lil_index'].value
    shape = lil_group.attrs['shape']
    lil = lil_matrix(np.empty(shape=shape), dtype=data.dtype)
    lil.data = data
    lil.rows = index
    if sparse_format == 'coo':
        return lil.tocoo()
    elif sparse_format == 'csr':
        return lil.tocsr()
    elif sparse_format == 'csc':
        return lil.tocsc()
    elif sparse_format == 'lil':
        return lil
    else:
        raise NotImplementedError('sparse_format only support coo, csr, csc, lil.')


def _get_mean_var(X):
    mean = X.mean(axis=0)
    if issparse(X):
        mean_sq = X.multiply(X).mean(axis=0)
        mean = mean.A1
        mean_sq = mean_sq.A1
    else:
        mean_sq = np.multiply(X, X).mean(axis=0)
    # enforece R convention (unbiased estimator) for variance
    var = (mean_sq - mean**2) * (X.shape[0]/(X.shape[0]-1))
    return mean, var


def _combine_mask(key, self_dict, obj_dict):
    if key in self_dict and key in obj_dict:
        new_row_mask = [all(l) for l in zip(self_dict[key], obj_dict[key])]
    elif key in self_dict:
        new_row_mask = self_dict
    elif key in obj_dict:
        new_row_mask = obj_dict
    else:
        new_row_mask = None
    return new_row_mask


def _slice_dict(attr_dict, select):
    # only select row
    # TODO change this function once change the dict into df
    df = pd.DataFrame(attr_dict)
    if isinstance(select, slice):
        return df.iloc[select, :].to_dict('series')
    elif isinstance(select, list):
        if isinstance(select[0], str):
            return df.loc[select, :].to_dict('series')
        else:
            return df.iloc[select, :].to_dict('series')
    else:
        raise NotImplementedError(f'Not support select type {type(select)}')


def combine_filter(s1, s2, axis=0):
    if axis in [0, 'row']:
        f1 = s1.rows['row_mask']
        f2 = s2.rows['row_mask']
        filter_bool = [all([x, y]) for x, y in zip(f1, f2)]
    else:
        f1 = s1.columns['col_mask']
        f2 = s2.columns['col_mask']
        filter_bool = list(np.concatenate([f1, f2]))
    return [idx for idx, filt in enumerate(filter_bool) if filt]

