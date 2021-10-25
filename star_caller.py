# -*- coding:utf-8 -*-
from collections import defaultdict
from argparse import ArgumentParser
from functools import reduce
import pandas as pd


DEGENERATE = {
    'R': 'A/G', 'Y': 'C/T', 'M': 'A/C', 'K': 'G/T', 'S': 'G/C', 'W': 'A/T'
}
COLUMN = 'Position at NC_000022.11 (Homo sapiens chromosome 22, GRCh38.p2)'


class StarCaller:
    def __init__(self, star_map):
        self.star_map = pd.read_excel(star_map, index_col=0)
        self.star_map = self.star_map[self.star_map['Structural Variation'].isna()]
        self.star_map.drop(columns=['Structural Variation'], inplace=True)
        self.star_map_dic = defaultdict(dict)
        for star, info in self.star_map.iterrows():
            self.star_map_dic[star] = info[~info.isna()].to_dict()

    def match(self, sample, result):
        columns = ['sampleID', 'star', COLUMN, 'g_hgvs_normalise', 'c_hgvs_normalise', 'star_allele', 'star_map', 'matched']
        matched_df = pd.DataFrame(columns=columns)
        result_df = pd.read_csv(result, sep='\t')
        result_df.index = result_df[COLUMN].tolist()
        result_dic = result_df.to_dict(orient='index')
        star_call = list()
        for star in self.star_map_dic:
            matched1_count, matched2_count = 0, 0
            for g in self.star_map_dic[star]:
                star_allele = str(self.star_map_dic[star][g])
                allele1, allele2 = result_dic[g]['star_allele'].split('/')
                if star_allele in DEGENERATE:
                    if allele1 in DEGENERATE[star_allele] and allele2 in DEGENERATE[star_allele]:
                        matched1_count += 1
                        matched2_count += 1
                        matched_result = 'yes/yes'
                    elif allele1 in DEGENERATE[star_allele] and allele2 not in DEGENERATE[star_allele]:
                        matched1_count += 1
                        matched_result = 'yes/no'
                    elif allele1 not in DEGENERATE[star_allele] and allele2 in DEGENERATE[star_allele]:
                        matched1_count += 1
                        matched_result = 'no/yes'
                    else:
                        matched_result = 'no/no'
                else:
                    if allele1 == star_allele and allele2 == star_allele:
                        matched1_count += 1
                        matched2_count += 1
                        matched_result = 'yes/yes'
                    elif allele1 == star_allele and allele2 != star_allele:
                        matched1_count += 1
                        matched_result = 'yes/no'
                    elif allele1 != star_allele and allele2 == star_allele:
                        matched1_count += 1
                        matched_result = 'no/yes'
                    else:
                        matched_result = 'no/no'
                matched_df = matched_df.append(pd.DataFrame(
                    [[sample, star, g, result_dic[g]['g_hgvs_normalise'], result_dic[g]['c_hgvs_normalise'],
                      result_dic[g]['star_allele'], star_allele, matched_result]], columns=columns), sort=False)
            if matched1_count == len(self.star_map_dic[star]) and matched2_count == len(self.star_map_dic[star]):
                star_call.append(star)
                star_call.append(star)
            elif matched1_count == len(self.star_map_dic[star]):
                star_call.append(star)
        star_call_dic = {sample: {'star': '/'.join(star_call)}}
        return star_call_dic, matched_df

    @staticmethod
    def star_adjust(star, cyp2d6):
        if cyp2d6 == '0':
            star_fixed = '*5/*5'
        elif cyp2d6 == '1':
            stars = set(star.split('/'))
            stars.add('*5')
            star_fixed = '/'.join(stars)
        else:
            star_fixed = star
        return star_fixed

    @classmethod
    def calling(cls, parsed_args):
        sc = cls(parsed_args.map)
        stars = list()
        df = pd.DataFrame()
        info = pd.read_csv(parsed_args.info, header=None, sep='\t')
        for i in range(info.shape[0]):
            sample, result = info.loc[i, 0], info.loc[i, 1]
            star_call_dic, matched_df = sc.match(sample, result)
            stars.append(star_call_dic)
            df = df.append(matched_df, sort=False)
        star_df = pd.DataFrame.from_dict(reduce(lambda x, y: dict(x, **y), stars), orient='index')
        star_df['sampleID'] = star_df.index.tolist()
        # 读取CYP2D6拷贝数
        try:
            cyp2d6_df = pd.read_csv(parsed_args.cyp2d6, sep='\t', dtype={'cyp2d6': str})
            star_df = pd.merge(star_df, cyp2d6_df, on=['sampleID'], how='left')
        except FileNotFoundError:
            star_df['cyp2d6'] = '.'
        star_df['star_fixed'] = star_df.apply(lambda x: sc.star_adjust(x['star'], x['cyp2d6']), axis=1)
        star_df.to_csv(f'{parsed_args.prefix}.cyp2d6_star.tsv', sep='\t', index=False, columns=['sampleID', 'star', 'star_fixed'])
        df.to_csv(f'{parsed_args.prefix}.cyp2d6_star_match.tsv', sep='\t', index=False)


def main():
    parser = ArgumentParser()
    parser.add_argument('-info', help='star allele parsed info', required=True)
    parser.add_argument('-map', help='CYP2D6_allele_definition_table map', required=True)
    parser.add_argument('-prefix', help='output file prefix', required=True)
    parser.add_argument('-cyp2d6', help='CYP2D6 copy number result file')
    parsed_args = parser.parse_args()
    StarCaller.calling(parsed_args)


if __name__ == '__main__':
    main()
