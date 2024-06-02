#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# @Time    : 2022/09/30
# @Author  : zhangsiwen
# @contact : zhang.siwen
# @ChangeLog
#     20220930, modified from /mnt/user/gaocx/pipeline/script/bin/bam_stat.py

import argparse
import json

# Get arguments
parser = argparse.ArgumentParser(description="Summary bam info from picard result.")
parser.add_argument("--fastp", "-f", help="fastp summary json", required=True)
parser.add_argument("--output", "-o", help="Out file, TSV format.", required=True)
args = parser.parse_args()


jsonfile = args.fastp
output = args.output

# read fastp summary json
fastp = json.load(open(jsonfile))

out_list = []
out_list.append("raw_reads\t%s" % fastp['summary']['before_filtering']['total_reads'])
out_list.append("raw_base\t%s" % fastp['summary']['before_filtering']['total_bases'])
out_list.append("clean_reads\t%s" % fastp['summary']['after_filtering']['total_reads'])
out_list.append("clean_bases\t%s" % fastp['summary']['after_filtering']['total_bases'])
ratio = 100 * float(fastp['summary']['after_filtering']['total_reads']) / float(fastp['summary']['before_filtering']['total_reads'])
out_list.append("clean_ratio\t{:.2f}%".format(ratio))
out_list.append("Q20_rate\t%s" % fastp['summary']['after_filtering']['q20_rate'])
out_list.append("Q30_rate\t%s" % fastp['summary']['after_filtering']['q30_rate'])
out_list.append("gc_content\t%s" % fastp['summary']['after_filtering']['gc_content'])

with open(output, "w") as f:
    f.writelines("\n".join(out_list) + "\n")
