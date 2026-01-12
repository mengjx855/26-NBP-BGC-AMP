#!/data/software/miniconda3/envs/jinxin/bin/python
# Jinxin Meng, mengjx855@163.com
# created date: 2025-09-29, 17:07:39
# modified date: 2025-11-12, 08:53:43

# 20251110 update: add multiprocessing.
# 20251112 update: select representative peptide not length-dependent, but MIC-value-dependent.

import re, sys, math, os, subprocess
from multiprocessing import Pool
from Bio import SeqIO
from collections import defaultdict

if len(sys.argv) != 4:
    sys.exit(f'Usage: {sys.argv[0]} [pep.fas] [avg.MIC] [out_prefix]')

TMP = f'tmp_{os.getpid()}'
SSW = '/data/software/SSWLibrary-1.2.5/src/ssw_test'
SIZE = 100
PROCESS = 10
RE_TARGET = re.compile(r'target_name: (\S+)')
RE_QUERY = re.compile(r'query_name: (\S+)')
RE_SCORE = re.compile(r'optimal_alignment_score: (\d+)')

if not os.path.isdir(TMP):
    os.makedirs(TMP)

# split fasta
fasta = sys.argv[1]
split_files = []
records = list(SeqIO.parse(fasta, 'fasta'))
width = len(str(math.ceil(len(records) / SIZE)))
for i in range(0, len(records), SIZE):
    file_num = i // SIZE + 1
    output_file = os.path.join(TMP, f'query.{file_num:0{width}d}')
    batch_records = records[i: i + SIZE]
    SeqIO.write(batch_records, output_file + '.fa', 'fasta')
    split_files.append(output_file) # split_files 是每个分块的前缀

# ssw_test
def run_ssw_command(args):
    target_fasta, query_prefix = args
    cmd = [SSW, '-p', target_fasta, f'{query_prefix}.fa']
    output_file = f'{query_prefix}.ssw'
    with open(output_file, 'w') as out_f, open(os.devnull, 'w') as err_f:
        subprocess.run(cmd, stdout=out_f, stderr=err_f, check=True) # check=True 如果返回码非0，自动抛出 CalledProcessError

args_list = [(fasta, split_file) for split_file in split_files]
with Pool(processes=PROCESS) as pool:
    pool.map(run_ssw_command, args_list)

# get scores
def recursive_dict(d):
    if isinstance(d, defaultdict):
        return {k: recursive_dict(v) for k, v in d.items()}
    return d

def run_get_scores(file_prefix):
    scores = defaultdict(lambda: defaultdict(lambda: 0))
    with open(file_prefix + '.ssw', 'r') as in_file:
        chunk = []
        for row in in_file:
            chunk.append(row.strip())
            if len(chunk) == 4:
                target = RE_TARGET.search(chunk[0]).group(1)
                query = RE_QUERY.search(chunk[1]).group(1)
                score = RE_SCORE.search(chunk[2]).group(1)
                scores[query][target] = int(score)
                chunk = []
    os.remove(file_prefix + '.ssw')
    os.remove(file_prefix + '.fa')
    return recursive_dict(scores)

args_list = [split_file for split_file in split_files]
with Pool(processes=PROCESS) as pool:
    results = pool.map(run_get_scores, args_list)
scores = defaultdict(dict) # 合并多进程的dict
for x in results:
    for k, v in x.items():
        scores[k] = v

# get_length
# length = {}
# seq_id = []
# for record in SeqIO.parse(fasta, 'fasta'):
#     length[record.id] = len(record.seq)
#     seq_id.append(record.id)

# MIC value
MIC = {}
in_file = open(sys.argv[2], 'r')
for row in in_file:
    parts = row.strip().split()
    MIC[parts[0]] = float(parts[1])

# calcu similarity and filter
seq_id = list(scores.keys())
out_file = open(sys.argv[3] + '.similarity', 'w')
similarity_score = defaultdict(lambda: defaultdict(lambda: 0))
unique = set(seq_id)
for x in range(len(seq_id)-1):
    for y in range(x + 1, len(seq_id)):
        x_ = seq_id[x]
        y_ = seq_id[y]
        score = scores[x_][y_] / math.sqrt(scores[x_][x_] * scores[y_][y_])
        similarity_score[x_][y_] = score
        out_file.write(f'{x_}\t{y_}\t{score}\n')
        if score > 0.7:
            if MIC[x_] <= MIC[y_]:
                unique.discard(y_)
            else:
                unique.discard(x_)
            # if length[x_] >= length[y_]:
            #     unique.discard(x_)
            # else:
            #     unique.discard(y_)
out_file.close()

out_file = open(sys.argv[3] + '.unique', 'w')
for record in SeqIO.parse(fasta, 'fasta'):
    if record.id in unique:
        out_file.write(f'{record.id}\t{record.seq}\n')
out_file.close()

os.rmdir(TMP)
