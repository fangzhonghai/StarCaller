# -*- coding:utf-8 -*-
from hgvs.assemblymapper import AssemblyMapper
from hgvs.dataproviders.uta import connect
from hgvs.normalizer import Normalizer
from argparse import ArgumentParser
from hgvs.parser import Parser
import pandas as pd
import pyfaidx
import re


COLUMN = 'Position at NC_000022.11 (Homo sapiens chromosome 22, GRCh38.p2)'


class Cyp2d6AlleleDefinition:
    def __init__(self, allele):
        self.allele = pd.read_excel(allele)
        self.allele.index = self.allele[COLUMN].tolist()

    @staticmethod
    def match(g, reference):
        if 'del' in g:
            pos, base = re.findall(r'g\.(\d+)del(\w+)', g)[0]
            position = int(pos) - 1
            base_front = str(reference.get_seq('chr22', position, position)).upper()
            ref_vcf = base_front + base
            alt_vcf = base_front
            ref_allele = base
            alt_allele = 'del' + base
            if len(base) == 1:
                g_hgvs = f'NC_000022.11:g.{pos}del{base}'
            else:
                g_hgvs = f'NC_000022.11:g.{pos}_{position+len(base)}del{base}'
        elif 'ins' in g:
            pos, base = re.findall(r'g\.(\d+).*ins(\w+)', g)[0]
            position = int(pos)
            base_front = str(reference.get_seq('chr22', position, position)).upper()
            ref_vcf = base_front
            alt_vcf = base_front + base
            ref_allele = 'del'
            alt_allele = 'ins' + base
            g_hgvs = f'NC_000022.11:{g}'
        else:
            pos, ref_vcf, alt_vcf = re.findall(r'g\.(\d+)(\w)>(\w)', g)[0]
            position = int(pos)
            ref_allele = ref_vcf
            alt_allele = alt_vcf
            g_hgvs = f'NC_000022.11:{g}'
        return str(position), ref_vcf, alt_vcf, ref_allele, alt_allele, g_hgvs

    @staticmethod
    def parsing(var_g, reference):
        var_g_list = var_g.split('; ')
        if len(var_g_list) == 1:
            position, ref_vcf, alt_vcf, ref_allele, alt_allele, g_hgvs = Cyp2d6AlleleDefinition.match(var_g_list[0], reference)
        else:
            parsed = list()
            for g in var_g_list:
                parsed.append((Cyp2d6AlleleDefinition.match(g, reference)))
            position, ref_vcf, alt_vcf, ref_allele, alt_allele, g_hgvs = (','.join(i) for i in zip(*parsed))
            position = position.split(',')[0]
            ref_vcf = ref_vcf.split(',')[0]
            ref_allele = ref_allele.split(',')[0]
        return position, ref_vcf, alt_vcf, ref_allele, alt_allele, g_hgvs

    @staticmethod
    def get_star_1(star_map_excel):
        star_map = pd.read_excel(star_map_excel, index_col=0)
        star_map = star_map[star_map['Structural Variation'].isna()].copy()
        star_map.drop(columns=['Structural Variation'], inplace=True)
        star_1 = star_map.loc['*1']
        star_1.name = 'ref_allele_star_1'
        return star_1

    @classmethod
    def process(cls, parsed_args):
        reference = pyfaidx.Fasta(parsed_args.reference)
        cad = cls(parsed_args.allele)
        cad.allele['parsing'] = cad.allele.apply(lambda x: ';'.join(cad.parsing(x[COLUMN], reference)), axis=1)
        df_split = cad.allele['parsing'].str.split(';', expand=True)
        df_split.columns = ['POS_origin', 'REF_origin', 'ALT_origin', 'ref_allele', 'alt_allele', 'g_hgvs']
        df_split = df_split.join(cad.get_star_1(parsed_args.map))
        df_out = cad.allele.join(df_split[['POS_origin', 'REF_origin', 'ALT_origin', 'ref_allele', 'ref_allele_star_1', 'alt_allele', 'g_hgvs']])
        df_out.drop(columns=['parsing'], inplace=True)
        df_out = normalise(df_out, reference)
        df_out['mut_type'] = 'snp'
        df_out.loc[df_out[COLUMN].str.contains('ins'), 'mut_type'] = 'ins'
        df_out.loc[df_out[COLUMN].str.contains('del'), 'mut_type'] = 'del'
        df_out.to_excel(parsed_args.out, index=False)


def annotator(annotation='GRCh38'):
    # 连接UTA, 创建Mapper、Parser、Normalizer
    hdp = connect()
    hp = Parser()
    # hn3 = Normalizer(hdp, shuffle_direction=3)
    hn5 = Normalizer(hdp, shuffle_direction=5)
    am = AssemblyMapper(hdp, assembly_name=annotation, alt_aln_method='splign', replace_reference=True)
    return hp, hn5, am


def normalise(df, reference):
    transcript = 'NM_000106.5'
    hp, hn5, am = annotator()
    df['g_hgvs_normalise'], df['c_hgvs_normalise'], df['POS_normalise'], df['REF_normalise'], df['ALT_normalise'] = '.', '.', '.', '.', '.'
    for i in df.index:
        g_hgvs_list = df.loc[i, 'g_hgvs'].split(',')
        normalised = list()
        for g in g_hgvs_list:
            g_parse = hp.parse_hgvs_variant(g)
            g_parse_normalise = hn5.normalize(g_parse)
            c_normalise = am.g_to_t(g_parse_normalise, transcript)
            if len(g_hgvs_list) == 1:
                position, ref_vcf, alt_vcf = coordinate_converter(g_parse_normalise, reference)
            else:
                normalised.append((str(g_parse_normalise), str(c_normalise), *coordinate_converter(g_parse_normalise, reference)))
        if normalised:
            g_parse_normalise, c_normalise, position, ref_vcf, alt_vcf = (','.join(i) for i in zip(*normalised))
            position = position.split(',')[0]
            ref_vcf = ref_vcf.split(',')[0]
        df.loc[i, 'g_hgvs_normalise'] = g_parse_normalise
        df.loc[i, 'c_hgvs_normalise'] = c_normalise
        df.loc[i, 'POS_normalise'] = position
        df.loc[i, 'REF_normalise'] = ref_vcf
        df.loc[i, 'ALT_normalise'] = alt_vcf
    return df


def coordinate_converter(g_parse, reference):
    if g_parse.posedit.edit.type == 'dup':
        position = g_parse.posedit.pos.start.base - 1
        base = str(reference.get_seq('chr22', position, position)).upper()
        ref_vcf = base
        alt_vcf = base + g_parse.posedit.edit.ref
    elif g_parse.posedit.edit.type == 'ins':
        position = g_parse.posedit.pos.start.base
        base = str(reference.get_seq('chr22', position, position)).upper()
        ref_vcf = base
        alt_vcf = base + g_parse.posedit.edit.alt
    elif g_parse.posedit.edit.type == 'del':
        position = g_parse.posedit.pos.start.base - 1
        base = str(reference.get_seq('chr22', position, position)).upper()
        ref_vcf = base + g_parse.posedit.edit.ref
        alt_vcf = base
    else:
        position = g_parse.posedit.pos.start.base
        ref_vcf = g_parse.posedit.edit.ref
        alt_vcf = g_parse.posedit.edit.alt
    return str(position), ref_vcf, alt_vcf


def main():
    parser = ArgumentParser()
    parser.add_argument('-reference', help='reference genome fasta file (hg38)', required=True)
    parser.add_argument('-allele', help='CYP2D6_allele_definition_table excel file', required=True)
    parser.add_argument('-map', help='CYP2D6_allele_definition_table map', required=True)
    parser.add_argument('-out', help='output excel file', required=True)
    parsed_args = parser.parse_args()
    Cyp2d6AlleleDefinition.process(parsed_args)


if __name__ == '__main__':
    main()
