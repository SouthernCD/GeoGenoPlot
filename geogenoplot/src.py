import matplotlib.pyplot as plt
import numpy as np
from scipy.stats import pearsonr, spearmanr, kruskal
from cyvcf2 import VCF
import pandas as pd
from yxmap import get_raster_score_many, scatter_map_ploter
import seaborn as sns
import warnings
warnings.filterwarnings('ignore')


def read_vcf_file(vcf_file, variant_id):
    """
    Read a VCF file and return a pandas.DataFrame with columns: 'CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO', 'FORMAT', 'SAMPLE'
    """
    vcf = VCF(vcf_file)
    records = []
    for record in vcf:
        if record.ID == variant_id:
            records.append(record)
    s_list = vcf.samples
    vcf.close()

    if records == []:
        raise ValueError("No variant with ID: %s" % variant_id)

    record = records[0]
    genotype = np.array(record.genotypes)

    var_df = pd.DataFrame(columns=['sample_id', 'genotype', 'phased'])
    for i, s in enumerate(s_list):
        if genotype[i][2]:
            genotype_str = '|'.join([str(x) for x in genotype[i][:2]])
        else:
            genotype_str = '/'.join([str(x) for x in genotype[i][:2]])
        var_df.loc[i] = [s, genotype_str, genotype[i][2]]

    return var_df, record.REF, record.ALT[0]


def add_location(var_df, coord_df):
    """
    Add 'latitude' and 'longitude' columns to var_df
    """
    for i, row in var_df.iterrows():
        var_df.loc[i, 'latitude'] = coord_df[coord_df['sample_id'] == row['sample_id']
                                             ]['latitude'].values[0] if row['sample_id'] in coord_df['sample_id'].values else np.nan
        var_df.loc[i, 'longitude'] = coord_df[coord_df['sample_id'] == row['sample_id']
                                              ]['longitude'].values[0] if row['sample_id'] in coord_df['sample_id'].values else np.nan

    return var_df


def add_raster_value(var_df, tif_file, factor=1):
    coord_dict = {}

    for i, row in var_df.iterrows():
        if not np.isnan(row['latitude']):
            coord_dict[row['sample_id']] = (row['latitude'], row['longitude'])

    value_dict = get_raster_score_many(coord_dict, tif_file)

    for i, row in var_df.iterrows():
        if var_df.loc[i, 'sample_id'] in value_dict.keys() and value_dict[var_df.loc[i, 'sample_id']] is not None:
            var_df.loc[i, 'value'] = value_dict[var_df.loc[i,
                                                           'sample_id']] * factor
        else:
            var_df.loc[i, 'value'] = np.nan

    return var_df


def allele_corr_plot(input_df, ref_allele=None, alt_allele=None, phased=False, trendline=True, y_label='Value', statistic='pearson', title='Correlation', save_file=None, ax=None):
    """
    input_df: pandas.DataFrame
    should have columns: genotype, and value columns
    genotype: genotype of the site (like '0|0', '0|1', '1|1'), 0 for reference allele, 1 for alternative allele
    value: value to be plotted as y-axis

    statistic = 'pearson', 'spearman' or 'kruskal'
    """
    input_df = input_df[input_df['value'].notnull()]
    input_df = input_df[input_df['genotype'].notnull()]

    if ax is None:
        fig, ax = plt.subplots(figsize=(4, 6))

    if phased:
        x_label_list = ['0|0', '0|1', '1|0', '1|1']
    else:
        x_label_list = ['0/0', '0/1', '1/1']
        ref_df = input_df.loc[input_df['genotype'].isin(['0|0', '0/0'])].copy()
        ref_df.loc[:, 'genotype'] = '0/0'
        hap_df = input_df.loc[input_df['genotype'].isin(
            ['0|1', '1|0', '1/0', '0/1'])].copy()
        hap_df.loc[:, 'genotype'] = '0/1'
        alt_df = input_df.loc[input_df['genotype'].isin(['1|1', '1/1'])].copy()
        alt_df.loc[:, 'genotype'] = '1/1'
        input_df = pd.concat([ref_df, hap_df, alt_df], ignore_index=True)

    sns.boxplot(x='genotype', y='value', data=input_df, showfliers=False,
                order=x_label_list, width=0.3, color='#BEDDFD', ax=ax)
    sns.stripplot(x='genotype', y='value', data=input_df, color='#258FF8',
                  jitter=0.2, order=x_label_list, size=3, alpha=0.5, ax=ax)

    # x = 0
    # x_label_list = []
    # x_list = []
    # y_list = []

    # df = input_df[input_df['genotype'].isin(['0|0', '0/0'])]
    # ax.scatter(np.full(len(df['value']), x), df['value'],
    #             color='#3B75AF', label='0|0' if phased else '0/0')
    # x_label_list.append('0|0' if phased else '0/0')
    # x_list.extend(list(np.full(len(df['value']), x)))
    # y_list.extend(list(df['value']))
    # x += 1

    # if phased:
    #     df = input_df[input_df['genotype'].isin(['0|1'])]
    #     ax.scatter(np.full(len(df['value']), x),
    #                 df['value'], color='#3B75AF', label='0|1')
    #     x_list.extend(list(np.full(len(df['value']), x)))
    #     y_list.extend(list(df['value']))
    #     x_label_list.append('0|1')
    #     x += 1

    #     df = input_df[input_df['genotype'].isin(['1|0'])]
    #     ax.scatter(np.full(len(df['value']), x),
    #                 df['value'], color='#3B75AF', label='1|0')
    #     x_list.extend(list(np.full(len(df['value']), x)))
    #     y_list.extend(list(df['value']))
    #     x_label_list.append('1|0')
    #     x += 1

    # else:
    #     df = input_df[input_df['genotype'].isin(['0|1', '1|0', '1/0', '0/1'])]
    #     ax.scatter(np.full(len(df['value']), x),
    #                 df['value'], color='#3B75AF', label='0/1')
    #     x_list.extend(list(np.full(len(df['value']), x)))
    #     y_list.extend(list(df['value']))
    #     x_label_list.append('0/1')
    #     x += 1

    # df = input_df[input_df['genotype'].isin(['1|1', '1/1'])]
    # ax.scatter(np.full(len(df['value']), x), df['value'],
    #             color='#3B75AF', label='1|1' if phased else '1/1')
    # x_list.extend(list(np.full(len(df['value']), x)))
    # y_list.extend(list(df['value']))
    # x_label_list.append('1|1' if phased else '1/1')
    # x += 1

    # ax.set_xticks(range(x))
    # ax.set_xticklabels(x_label_list)

    ax.set_xlim(-0.5, len(x_label_list)-0.5)

    if ref_allele is not None and alt_allele is not None:
        x_label = "Genotype (Ref: {0}, Alt: {1})".format(
            ref_allele, alt_allele)
    else:
        x_label = "Genotype"

    ax.set_xlabel(x_label)
    ax.set_ylabel(y_label)

    if not phased:
        # 拟合数据并绘制趋势线
        if trendline:
            x_list = input_df['genotype']
            x_list.replace('0/0', 0, inplace=True)
            x_list.replace('0/1', 1, inplace=True)
            x_list.replace('1/1', 2, inplace=True)
            y_list = input_df['value']

            try:
                z = np.polyfit(x_list, y_list, 1)
                p = np.poly1d(z)
                ax.plot(x_list, p(x_list), "-", c='#E87D85')
            except:
                pass

        if statistic == 'pearson':
            p = pearsonr(x_list, y_list)
            subtitile = "Pearson correlation: %.4e, p-value: %.4e" % (
                p[0], p[1])
        elif statistic == 'spearman':
            s = spearmanr(x_list, y_list)
            subtitile = "Spearman correlation: %.4e, p-value: %.4e" % (
                s.correlation, s.pvalue)
        elif statistic == 'kruskal':
            k = kruskal(x_list, y_list)
            subtitile = "Kruskal statistic: %.4e, p-value: %.4e" % (
                k.statistic, k.pvalue)

        ax.text(0.5, 1.12, title, transform=ax.transAxes,
                ha="center", va="center", fontsize=18)
        ax.text(0.5, 1.05, subtitile, transform=ax.transAxes,
                ha="center", va="center", fontsize=12)
    else:
        ax.text(0.5, 1.05, title, transform=ax.transAxes,
                ha="center", va="center", fontsize=18)

    if ax is None:
        plt.show()

    if save_file is not None:
        fig.savefig(save_file, format='pdf', facecolor='none',
                    edgecolor='none', bbox_inches='tight')
    return ax


def plot_main(args):

    var_df, ref, alt = read_vcf_file(args.vcf_file, args.var_id)
    coord_df = pd.read_csv(args.georef_table)
    var_df = add_location(var_df, coord_df)
    var_df = add_raster_value(var_df, args.raster_file)

    if args.statistic == 'all':
        allele_corr_plot(var_df, ref_allele=ref, alt_allele=alt, phased=args.phased, trendline=not args.no_trendline,
                         y_label=args.env_factor_name, statistic='pearson', title='Correlation', save_file="%s.correlation.peason.pdf" % args.var_id)
        allele_corr_plot(var_df, ref_allele=ref, alt_allele=alt, phased=args.phased, trendline=not args.no_trendline,
                         y_label=args.env_factor_name, statistic='spearman', title='Correlation', save_file="%s.correlation.spearman.pdf" % args.var_id)
        allele_corr_plot(var_df, ref_allele=ref, alt_allele=alt, phased=args.phased, trendline=not args.no_trendline,
                         y_label=args.env_factor_name, statistic='kruskal', title='Correlation', save_file="%s.correlation.kruskal.pdf" % args.var_id)
    else:
        allele_corr_plot(var_df, ref_allele=ref, alt_allele=alt, phased=args.phased, trendline=not args.no_trendline,
                         y_label=args.env_factor_name, statistic=args.statistic, title='Correlation', save_file="%s.correlation.%s.pdf" % (args.var_id, args.statistic))

    var_df.to_csv('%s.csv' % args.var_id, index=False)

    data = var_df[['latitude', 'longitude', 'genotype', 'value']]
    data.columns = ['latitude', 'longitude', 'type', 'value']
    data

    defaut_longitude_limit = [int(i) for i in args.lon.split(',')]
    defaut_latitude_limit = [int(i) for i in args.lat.split(',')]

    fig, ax = plt.subplots(figsize=(15, 8))

    scatter_map_ploter(data,
                       title=args.var_id,
                       value_label=args.env_factor_name,
                       longitude_limit=defaut_longitude_limit,
                       latitude_limit=defaut_latitude_limit,
                       type_plot_order=['0|0', '0|1', '1|0', '1|1'],
                       ax=ax,
                       s=20)

    fig.savefig('%s.map.pdf' % args.var_id, format='pdf', facecolor='none',
                edgecolor='none', bbox_inches='tight')


if __name__ == '__main__':
    pass
