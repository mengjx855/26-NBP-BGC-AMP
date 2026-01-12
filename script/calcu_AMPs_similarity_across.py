#!/data/software/miniconda3/envs/jinxin/bin/python
# Jinxin Meng, mengjx855@163.com
# created date: 2025-09-29, 17:07:39
# modified date: 2025-11-10, 20:04:36

import re, sys, math, os, subprocess
from multiprocessing import Pool
from Bio import SeqIO
from collections import defaultdict

if len(sys.argv) != 4:
    sys.exit(f'Usage: {sys.argv[0]} [pep.fas] [ref.fas] [out_prefix]')

TMP = f'tmp_{os.getpid()}'
SSW = '/data/software/SSWLibrary-1.2.5/src/ssw_test'
PROCESS = 10
SIZE = 100
RE_TARGET = re.compile(r'target_name: (\S+)')
RE_QUERY = re.compile(r'query_name: (\S+)')
RE_SCORE = re.compile(r'optimal_alignment_score: (\d+)')

if not os.path.isdir(TMP):
    os.makedirs(TMP)

# cross, split query fasta
query = sys.argv[1]
target = sys.argv[2]
split_files = []
records = list(SeqIO.parse(query, 'fasta'))
width = len(str(math.ceil(len(records) / SIZE)))
for i in range(0, len(records), SIZE):
    file_num = i // SIZE + 1
    output_file = os.path.join(TMP, f'query.{file_num:0{width}d}')
    batch_records = records[i: i + SIZE]
    SeqIO.write(batch_records, output_file + '.fa', 'fasta')
    split_files.append(output_file) # split_files 是每个分块的前缀

# cross, ssw_test with target
def run_ssw_command(args):
    target_fasta, query_prefix = args
    cmd = [SSW, '-p', target_fasta, f'{query_prefix}.fa']
    output_file = f'{query_prefix}.ssw'
    with open(output_file, 'w') as out_f, open(os.devnull, 'w') as err_f:
        subprocess.run(cmd, stdout=out_f, stderr=err_f, check=True) # check=True 如果返回码非0，自动抛出 CalledProcessError
args_list = [(target, split_file) for split_file in split_files]
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
across_scores = defaultdict(dict) # 合并多进程的dict
for x in results:
    for k, v in x.items():
        across_scores[k] = v

# self ssw and get scores
def run_ssw_self(record):
    output_file = os.path.join(TMP, f'.self.{record.id}.fa')
    SeqIO.write(record, output_file, 'fasta')
    result = subprocess.run([SSW, '-p', output_file, output_file], stdout=subprocess.PIPE, stderr=subprocess.DEVNULL, text=True, check=True)
    output = result.stdout.replace('\n',' ')
    score = RE_SCORE.search(output).group(1)
    os.remove(output_file)
    return (record.id, score)

self_scores = defaultdict(lambda: 0)
records = list(SeqIO.parse(query, 'fasta')) # query_self, -> self_scores
with Pool(processes=PROCESS) as pool:
    results = pool.map(run_ssw_self, records)
for x in results:
    self_scores[x[0]] = int(x[1])
query_id = []
for record in records:
    query_id.append(record.id)

records = list(SeqIO.parse(target, 'fasta')) # target_self, -> self_scores
with Pool(processes=PROCESS) as pool:
    results = pool.map(run_ssw_self, records)
for x in results:
    self_scores[x[0]] = int(x[1])
target_id = []
for record in records:
    target_id.append(record.id)

# calcu similarity and filter
out_file = open(sys.argv[3] + '.similarity', 'w')
unique = set(query_id)
flag = False
for x in query_id:
    for y in target_id:
        score = across_scores[x][y] / math.sqrt(self_scores[x] * self_scores[y])
        out_file.write(f'{x}\t{y}\t{score}\n')
        if score > 0.7:
            flag = True
    if flag:
        unique.discard(x)
    flag = False
out_file.close()

out_file = open(sys.argv[3] + '.unique', 'w')
for record in SeqIO.parse(query, 'fasta'):
    if record.id in unique:
        out_file.write(f'{record.id}\t{record.seq}\n')
out_file.close()

os.rmdir(TMP)
