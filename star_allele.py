# -*- coding:utf-8 -*-
from collections import defaultdict, namedtuple
from argparse import ArgumentParser
from pysam import AlignmentFile
import pandas as pd
import numpy as np
import tabix
import re
import os


VariantGenotype = namedtuple('VariantGenotype', ['ref_calling', 'alt_calling', 'QD', 'GT', 'DP', 'Ratio', 'genotype',
                                                 'GT_corrected', 'genotype_corrected', 'star_allele'])
VariantRecord = namedtuple('VariantRecord', ['ref_calling', 'alt_calling', 'QD', 'GT', 'DP', 'Ratio'])
ROOT = os.path.abspath(os.path.dirname(__file__))
CONFIG = os.path.join(ROOT, 'etc')


class StarAllele:
    def __init__(self, bam, vcf):
        self.bam = AlignmentFile(bam, 'rb')
        self.vcf = tabix.open(vcf)
        self.GenotypeDict = defaultdict(dict)

    @staticmethod
    def position_depth(bam, chromosome, position):
        depth = np.sum(bam.count_coverage(contig=chromosome, start=position - 1, stop=position))
        return depth

    def parsing(self, chromosome, position):
        try:
            vcf_record = next(self.vcf.query(chromosome, position - 1, position))
        except (StopIteration, tabix.TabixError):
            gt = '0/0'
            dp = StarAllele.position_depth(self.bam, chromosome, position)
            return VariantRecord('.', '.', '.', gt, dp, '.')
        ref, alt = vcf_record[3], vcf_record[4]
        try:
            qd = re.findall('QD=(.*?);', vcf_record[7])[0]
        except IndexError:
            qd = '.'
        info = vcf_record[9].split(':')
        gt = info[0].replace('|', '/')
        ad, dp = info[1], int(info[2])
        # ratio = int(ad.split(',')[-1]) / dp
        ratio = int(ad.split(',')[-1]) / np.sum([int(i) for i in ad.split(',')])
        return VariantRecord(ref, alt, qd, gt, dp, ratio)

    @staticmethod
    def genotyping(drug_record, parsed_record):
        # genotype_corrected,gt_corrected,star_allele,mut_type设置初始值
        genotype_corrected = 'Fail'
        gt_corrected = parsed_record.GT
        star_allele = f'{drug_record.ref_allele}/{drug_record.ref_allele}'
        mut_type = judge_mut_type(parsed_record.ref_calling, parsed_record.alt_calling)
        # 分情况讨论,目前gt_corrected都会根据ratio进行矫正;
        # drug_record.corrected控制genotype_corrected和star_allele是否为矫正后的;
        # drug_record.depth_threshold控制star_allele是否为'./.',genotype_corrected是否为'Fail'
        # 0/0
        if parsed_record.GT == '0/0':
            genotype = f'{drug_record.REF_normalise}/{drug_record.REF_normalise}'
            # 判断是否通过深度阈值进行 基因/star 型别分析,未赋值的变量为初始默认值
            if parsed_record.DP >= drug_record.depth_threshold:
                genotype_corrected = genotype
            else:
                star_allele = './.'
        # 0/1
        elif parsed_record.GT in ['0/1', '1/0']:
            genotype = f'{parsed_record.ref_calling}/{parsed_record.alt_calling}'
            # 判断是否通过深度阈值进行 基因/star 型别分析
            if parsed_record.DP >= drug_record.depth_threshold:
                if parsed_record.Ratio > drug_record.het_hom_ratio_threshold:
                    genotype_corrected = f'{parsed_record.alt_calling}/{parsed_record.alt_calling}'
                    gt_corrected = '1/1'
                elif parsed_record.Ratio < drug_record.ref_het_ratio_threshold:
                    genotype_corrected = f'{drug_record.REF_normalise}/{drug_record.REF_normalise}'
                    gt_corrected = '0/0'
                else:
                    genotype_corrected = genotype
                # 根据检测到的变异类型和是否使用矫正的基因型进行分情况讨论,不同变异类型star型别写法不同
                if mut_type == 'del':
                    star_allele = f'{drug_record.ref_allele}/del{parsed_record.ref_calling[1:]}'
                    if drug_record.corrected == 'yes':
                        if gt_corrected == '1/1':
                            star_allele = f'del{parsed_record.ref_calling[1:]}/del{parsed_record.ref_calling[1:]}'
                        elif gt_corrected == '0/0':
                            star_allele = f'{drug_record.ref_allele}/{drug_record.ref_allele}'
                elif mut_type == 'ins':
                    star_allele = f'del/ins{parsed_record.alt_calling[1:]}'
                    if drug_record.corrected == 'yes':
                        if gt_corrected == '1/1':
                            star_allele = f'ins{parsed_record.alt_calling[1:]}/ins{parsed_record.alt_calling[1:]}'
                        elif gt_corrected == '0/0':
                            star_allele = f'{drug_record.ref_allele}/{drug_record.ref_allele}'
                else:
                    star_allele = f'{parsed_record.ref_calling}/{parsed_record.alt_calling}'
                    if drug_record.corrected == 'yes':
                        if gt_corrected == '1/1':
                            star_allele = f'{parsed_record.alt_calling}/{parsed_record.alt_calling}'
                        elif gt_corrected == '0/0':
                            star_allele = f'{drug_record.ref_allele}/{drug_record.ref_allele}'
            else:
                star_allele = './.'
        # 1/1
        elif parsed_record.GT == '1/1':
            genotype = f'{parsed_record.alt_calling}/{parsed_record.alt_calling}'
            # 判断是否通过深度阈值进行 基因/star 型别分析
            if parsed_record.DP >= drug_record.depth_threshold:
                if parsed_record.Ratio <= drug_record.het_hom_ratio_threshold:
                    genotype_corrected = f'{parsed_record.ref_calling}/{parsed_record.alt_calling}'
                    gt_corrected = '0/1'
                else:
                    genotype_corrected = genotype
                # 根据检测到的变异类型和是否使用矫正的基因型进行分情况讨论,不同变异类型star型别写法不同
                if mut_type == 'del':
                    star_allele = f'del{parsed_record.ref_calling[1:]}/del{parsed_record.ref_calling[1:]}'
                    if drug_record.corrected == 'yes':
                        if gt_corrected == '0/1':
                            star_allele = f'{drug_record.ref_allele}/del{parsed_record.ref_calling[1:]}'
                elif mut_type == 'ins':
                    star_allele = f'ins{parsed_record.alt_calling[1:]}/ins{parsed_record.alt_calling[1:]}'
                    if drug_record.corrected == 'yes':
                        if gt_corrected == '0/1':
                            star_allele = f'del/ins{parsed_record.alt_calling[1:]}'
                else:
                    star_allele = f'{parsed_record.alt_calling}/{parsed_record.alt_calling}'
                    if drug_record.corrected == 'yes':
                        if gt_corrected == '0/1':
                            star_allele = f'{parsed_record.ref_calling}/{parsed_record.alt_calling}'
            else:
                star_allele = './.'
        # 1/2
        else:
            alt_1, alt_2 = (alt for alt in parsed_record.alt_calling.split(','))
            genotype = f'{alt_1}/{alt_2}'
            ratio_1, ratio_2 = 1 - parsed_record.Ratio, parsed_record.Ratio
            # 判断是否通过深度阈值进行 基因/star 型别分析
            if parsed_record.DP >= drug_record.depth_threshold:
                if ratio_1 > drug_record.het_hom_ratio_threshold:
                    genotype_corrected = f'{alt_1}/{alt_1}'
                    gt_corrected = '1/1'
                elif ratio_2 > drug_record.het_hom_ratio_threshold:
                    genotype_corrected = f'{alt_2}/{alt_2}'
                    gt_corrected = '1/1'
                else:
                    genotype_corrected = genotype
                # 根据检测到的变异类型和是否使用矫正的基因型进行分情况讨论,不同变异类型star型别写法不同
                mut_type_1 = judge_mut_type(parsed_record.ref_calling, alt_1)
                mut_type_2 = judge_mut_type(parsed_record.ref_calling, alt_2)
                if mut_type_1 == 'del' or mut_type_2 == 'del':
                    star_allele = f'del{parsed_record.ref_calling[len(ratio_1):]}/del{parsed_record.ref_calling[len(ratio_2):]}'
                    if drug_record.corrected == 'yes':
                        if ratio_1 > drug_record.het_hom_ratio_threshold:
                            star_allele = f'del{parsed_record.ref_calling[len(ratio_1):]}/del{parsed_record.ref_calling[len(ratio_1):]}'
                        elif ratio_2 > drug_record.het_hom_ratio_threshold:
                            star_allele = f'del{parsed_record.ref_calling[len(ratio_2):]}/del{parsed_record.ref_calling[len(ratio_2):]}'
                elif mut_type_1 == 'ins' or mut_type_2 == 'ins':
                    star_allele = f'ins{alt_1[1:]}/ins{alt_2[1:]}'
                    if drug_record.corrected == 'yes':
                        if ratio_1 > drug_record.het_hom_ratio_threshold:
                            star_allele = f'ins{alt_1[1:]}/ins{alt_1[1:]}'
                        elif ratio_2 > drug_record.het_hom_ratio_threshold:
                            star_allele = f'ins{alt_2[1:]}/ins{alt_2[1:]}'
                else:
                    star_allele = f'{alt_1}/{alt_2}'
                    if drug_record.corrected == 'yes':
                        if ratio_1 > drug_record.het_hom_ratio_threshold:
                            star_allele = f'{alt_1}/{alt_1}'
                        elif ratio_2 > drug_record.het_hom_ratio_threshold:
                            star_allele = f'{alt_2}/{alt_2}'
            else:
                star_allele = './.'
        if drug_record.corrected == 'no' and genotype_corrected != 'Fail':
            genotype_corrected = genotype
        return VariantGenotype(*parsed_record, genotype, gt_corrected, genotype_corrected, star_allele)

    @classmethod
    def calling(cls, parsed_args):
        sa = cls(parsed_args.bam, parsed_args.vcf)
        column = 'Position at NC_000022.11 (Homo sapiens chromosome 22, GRCh38.p2)'
        drug = pd.read_excel(os.path.join(CONFIG, parsed_args.variant + '.xlsx'))
        drug.index = drug[column].tolist()
        for drug_record in drug.itertuples():
            parsed_record = sa.parsing('chr22', drug_record.POS_normalise)
            sa.GenotypeDict[drug_record.Index] = sa.genotyping(drug_record, parsed_record)._asdict()
        sa.bam.close()
        df = pd.DataFrame.from_dict(sa.GenotypeDict, orient='index')
        df = drug.join(df)
        df.to_csv(parsed_args.out, sep='\t', index=False)


def judge_mut_type(ref_calling, alt_calling):
    if len(ref_calling) > 1 and len(alt_calling) == 1:
        mut_type = 'del'
    elif len(ref_calling) == 1 and len(alt_calling) > 1:
        mut_type = 'ins'
    elif ref_calling == '.' or alt_calling == '.':
        mut_type = 'ref'
    elif len(ref_calling) == 1 and len(alt_calling) == 1:
        mut_type = 'snp'
    else:
        mut_type = 'delins'
    return mut_type


def main():
    parser = ArgumentParser()
    parser.add_argument('-bam', help='bam file', required=True)
    parser.add_argument('-vcf', help='vcf file', required=True)
    parser.add_argument('-out', help='output file', required=True)
    parser.add_argument('-variant', help='etc excel file prefix', default='variant')
    parsed_args = parser.parse_args()
    StarAllele.calling(parsed_args)


if __name__ == '__main__':
    main()
