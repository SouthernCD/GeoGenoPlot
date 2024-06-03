import argparse
from geogenoplot import plot_main

class Job(object):
    def __init__(self):
        pass

    def run_arg_parser(self):
        # argument parse
        parser = argparse.ArgumentParser(
            prog='geogenoplot',
        )

        # argparse for distance
        parser.add_argument('vcf_file', type=str,
                            help='input vcf file')
        parser.add_argument('raster_file', type=str,
                            help='input raster file')
        parser.add_argument('georef_table', type=str,
                            help='input georeference table, should be a csv file with columns: sample_id, latitude, longitude')
        parser.add_argument('var_id', type=str,
                            help='variable id, should be match ID column in the vcf file')
        parser.add_argument('-p', '--phased', action='store_true', default=False, help='phased or not, default is not phased')
        parser.add_argument('-n', '--env_factor_name', type=str,
                            help='environmental factor name', default='value')
        parser.add_argument('-s', '--statistic', type=str,
                            default='all', help='correlation statistic, choices are pearson, spearman, kruskal, all. Default is all.',
                            choices=['pearson', 'spearman', 'kruskal', 'all'])
        parser.add_argument('-t', '--no_trendline', action='store_true', default=False, help='no trendline, default is to show trendline')
        parser.add_argument('--lat', type=str, default='-37,51', help='latitude range, default is -37,51')
        parser.add_argument('--lon', type=str, default='-20,95', help='longitude range, default is -20,95')

        self.arg_parser = parser

        self.args = parser.parse_args()


    def run(self):
        self.run_arg_parser()

        if self.args.vcf_file and self.args.raster_file and self.args.georef_table and self.args.var_id:
            plot_main(self.args)
        else:
            self.arg_parser.print_help()

def main():
    job = Job()
    job.run()


if __name__ == '__main__':
    main()