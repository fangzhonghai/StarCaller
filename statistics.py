# -*- coding:utf-8 -*-
from collections import defaultdict
from argparse import ArgumentParser
import pandas as pd


def statistics(df):
    stat = defaultdict(dict)
    df = df[df['Structural Variation'].isna()].copy()
    df.drop(columns=['Structural Variation'], inplace=True)
    for star, info in df.iterrows():
        stat[star]['specific_variant'] = '|'.join(info[~info.isna()].index)
        stat[star]['specific_variant_count'] = sum(~info.isna())
        stat[star]['degenerate_variant'] = '|'.join(info[info.isin(['R', 'Y', 'M', 'K', 'S', 'W'])].index)
        stat[star]['degenerate_variant_count'] = sum(info.isin(['R', 'Y', 'M', 'K', 'S', 'W']))
    out = pd.DataFrame.from_dict(stat, orient='index')
    return out


def main():
    parser = ArgumentParser()
    parser.add_argument('-map', help='CYP2D6_allele_definition_table map', required=True)
    parser.add_argument('-out', help='output file', required=True)
    parsed_args = parser.parse_args()
    df = pd.read_excel(parsed_args.map, index_col=0)
    out = statistics(df)
    out.to_excel(parsed_args.out)


if __name__ == '__main__':
    main()
