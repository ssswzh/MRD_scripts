#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# @Author  : weibin & zhangsiwen
# @ChangeLog
#     20220522, first version

import datetime
import sys
import argparse
from openpyxl import load_workbook


def GetArgs():
    parser = argparse.ArgumentParser(description='fill sites to order template', formatter_class=argparse.RawTextHelpFormatter)
    parser._action_groups.pop()
    # required arguments
    required = parser.add_argument_group('Required arguments')
    required.add_argument('--tsv', help='sites tsv file, columns are: Hugo_Symbol\tChromosome\tStart_Position\tEnd_Position\tHGVSc\tHGVSp', action='store', dest='tsv', required=True)
    required.add_argument('--id', help='sample id to be displayed in order worksheet', action='store', dest='id', required=True)
    required.add_argument('--out', help='out file name.xlsx', action='store', dest='out', required=True)
    required.add_argument('--temp', help='pyclone input file', action='store', dest='temp', required=True)
    args = parser.parse_args()
    return args


def now():
    return datetime.datetime.now().strftime("%Y/%m/%d")


def generate_file(template_file, in_file, sampleid, out_file):
    workbook = load_workbook(template_file)
    worksheet = workbook['CUSTOM']
    worksheet['C4'] = now()
    worksheet['C5'] = sampleid

    with open(in_file, "r") as f:
        row_idx = 5
        for line in f:
            data = line.strip("\n").split("\t")

            worksheet["G%s" % row_idx] = data[0]
            worksheet["H%s" % row_idx] = data[1]
            worksheet["I%s" % row_idx] = data[2]
            worksheet["J%s" % row_idx] = data[3]
            worksheet["K%s" % row_idx] = data[4]
            worksheet["L%s" % row_idx] = data[5]

            row_idx += 1

    workbook.save(out_file)


def main():
    args = GetArgs()
    generate_file(args.temp, args.tsv, args.id, args.out)


if __name__ == '__main__':
    main()
