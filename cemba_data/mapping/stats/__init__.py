import pathlib
import subprocess

import pandas as pd
from papermill import execute_notebook, PapermillExecutionError

from .m3c import m3c_mapping_stats, m3c_additional_cols
from .mc import mc_mapping_stats, mc_additional_cols
from .mct import mct_mapping_stats, mct_additional_cols
from .plate_info import get_plate_info
from ..pipelines import PACKAGE_DIR
from ...utilities import get_configuration


def mapping_stats(output_dir):
    """This is UID level mapping summary, the config file is in parent dir"""
    output_dir = pathlib.Path(output_dir).absolute()
    config = get_configuration(output_dir.parent / 'mapping_config.ini')
    mode = config['mode']

    if mode == 'mc':
        final_df = mc_mapping_stats(output_dir, config)
    elif mode == 'mct':
        final_df = mct_mapping_stats(output_dir, config)
    elif mode == 'm3c':
        final_df = m3c_mapping_stats(output_dir, config)
    else:
        raise ValueError

    # plate info, which is tech independent.
    _plate_info = get_plate_info(final_df.index, barcode_version=config['barcode_version'])
    final_df = pd.concat([_plate_info, final_df], axis=1)

    # save
    final_df.to_csv(output_dir / 'MappingSummary.csv.gz')
    return


def final_summary(output_dir, cleanup=True):
    output_dir = pathlib.Path(output_dir).absolute()
    mode = get_configuration(output_dir / 'mapping_config.ini')['mode']
    path_to_remove = []

    # Before running summary,
    # first make sure all the UID dir having Snakefile also has mapping summary (means successful)
    snakefile_list = list(output_dir.glob('*/Snakefile'))
    summary_paths = []
    for path in snakefile_list:
        uid_dir = path.parent
        summary_path = uid_dir / 'MappingSummary.csv.gz'
        if summary_path.exists():
            summary_paths.append(summary_path)
        else:
            raise FileNotFoundError(f'All sub dir should be successfully mapped before generating final summary. \n'
                                    f'However, MappingSummary.csv.gz does not exist in sub dir {uid_dir}. \n'
                                    f'This is the final target file of snakefile in {path}. \n'
                                    f'Run the following command to see what is missing: \n'
                                    f'snakemake -d {uid_dir} --snakefile {path} -np')

    # aggregate mapping summaries
    total_mapping_summary = pd.concat([pd.read_csv(path, index_col=0)
                                       for path in summary_paths])
    total_mapping_summary_path = output_dir / 'stats/MappingSummary.csv.gz'

    # if this is mct, aggregate all the gene counts
    if mode == 'mct':
        from ..stats.mct import aggregate_feature_counts
        aggregate_feature_counts(output_dir)

    # add additional columns based on some calculation
    if mode == 'mc':
        total_mapping_summary = mc_additional_cols(total_mapping_summary)
    elif mode == 'mct':
        total_mapping_summary = mct_additional_cols(total_mapping_summary, output_dir=output_dir)
    elif mode == 'm3c':
        total_mapping_summary = m3c_additional_cols(total_mapping_summary)
    else:
        raise

    # save total mapping summary
    total_mapping_summary.to_csv(total_mapping_summary_path)

    # add .snakemake files to deletion
    snakemake_hiding_dirs = list(output_dir.glob('*/.snakemake'))
    path_to_remove += snakemake_hiding_dirs

    # run summary notebook
    nb_path = output_dir / 'stats/MappingSummary.ipynb'
    try:
        mode = get_configuration(output_dir / 'mapping_config.ini')['mode']
        template_notebook = PACKAGE_DIR / f'files/mapping_summary_template/{mode}_template.ipynb'
        print('Executing summary plotting notebook...')
        execute_notebook(
            input_path=str(template_notebook),
            output_path=str(nb_path),
            parameters=dict(output_dir=str(output_dir))
        )
        print('Summary notebook successfully executed. Exporting HTML...')
        subprocess.run(['jupyter', 'nbconvert', '--to', 'html', str(nb_path)])
        print(f'See the summary plots here: {str(nb_path)[:-5]}html')
        print(f'Or customize the summary plots here: {nb_path}')
    except PapermillExecutionError:
        print(f'Ops, summary plotting notebook got some error, check the information in {nb_path}')
        cleanup = False

    # delete
    if cleanup:
        print('Clean up snakemake log...')
        for path in path_to_remove:
            subprocess.run(['rm', '-rf', str(path)], check=True)
    return