import subprocess
import tempfile
import os
import glob
import time
import json
import sys
import re
import shutil
from io import StringIO
from Bio import PDB, SeqIO

def parse_fasta_block(lines):
    """从行列表（一个块）中解析FASTA记录，返回序列列表"""
    fasta_str = '\n'.join(lines)
    handle = StringIO(fasta_str)
    sequences = []
    for record in SeqIO.parse(handle, 'fasta'):
        sequences.append(str(record.seq))
    return sequences

def parse_batch_file(file_path):
    """
    解析批量输入文件，返回列表，每个元素为 (identifier, seq1, seq2)
    文件格式：
        ID1
        >chainA
        seqA...
        >chainB
        seqB...

        ID2
        >chainA
        ...
    """
    with open(file_path, 'r') as f:
        content = f.read()
    blocks = re.split(r'\n\s*\n', content.strip())
    cases = []
    for block in blocks:
        lines = block.strip().split('\n')
        if not lines:
            continue
        identifier = lines[0].strip()
        fasta_lines = lines[1:]
        if not fasta_lines:
            print(f"警告：块 {identifier} 没有序列数据，跳过")
            continue
        seqs = parse_fasta_block(fasta_lines)
        if len(seqs) < 2:
            print(f"警告：块 {identifier} 序列不足2条（实际{len(seqs)}），跳过")
            continue
        cases.append((identifier, seqs[0], seqs[1]))
    return cases

def predict_interchain_disulfide(seq1, seq2, chain1_id='A', chain2_id='B', num_recycle=3, distance_threshold=2.5, save_pdb_path=None):
    """核心预测函数，使用MSA模式（mmseqs2_uniref），可选择性保存PDB文件"""
    with tempfile.TemporaryDirectory() as tmpdir:
        fasta_path = os.path.join(tmpdir, 'input.fasta')
        with open(fasta_path, 'w') as f:
            f.write(f">{chain1_id}:{chain2_id}\n{seq1}:{seq2}\n")

        wsl_fasta = fasta_path.replace('C:', '/mnt/c').replace('\\', '/')
        wsl_outdir = tmpdir.replace('C:', '/mnt/c').replace('\\', '/')

        # 请根据你的WSL用户名修改 -u 后面的参数（默认为 jmwang24）
        cmd = (f'wsl -u jmwang24 bash -c "source /home/jmwang24/miniconda3/etc/profile.d/conda.sh && '
               f'conda activate colabfold && '
               f'colabfold_batch {wsl_fasta} {wsl_outdir} --num-models 1 --model-type alphafold2_multimer_v3 --num-recycle {num_recycle} --msa-mode mmseqs2_uniref"')

        print(f"  正在运行 ColabFold 预测...")
        sys.stdout.flush()
        start = time.time()
        process = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
                                   text=True, bufsize=1, universal_newlines=True, encoding='utf-8', errors='ignore')
        for line in process.stdout:
            print(line, end='')
            sys.stdout.flush()
        process.wait()
        if process.returncode != 0:
            print("  预测失败")
            return {"disulfides": [], "pdb_file": None, "sequences": {"chain1": seq1, "chain2": seq2}}

        print(f"  预测完成，耗时 {time.time()-start:.1f} 秒")

        # 查找PDB文件
        pdb_files = glob.glob(os.path.join(tmpdir, '*.pdb'))
        if not pdb_files:
            print("  未找到PDB文件")
            return {"disulfides": [], "pdb_file": None, "sequences": {"chain1": seq1, "chain2": seq2}}
        ranked_pdbs = [f for f in pdb_files if 'rank_1' in f]
        pdb_file = ranked_pdbs[0] if ranked_pdbs else pdb_files[0]

        # 保存PDB到指定路径（如果提供）
        if save_pdb_path:
            os.makedirs(os.path.dirname(save_pdb_path), exist_ok=True)
            shutil.copy(pdb_file, save_pdb_path)
            print(f"  PDB已保存至 {save_pdb_path}")

        # 解析PDB，计算距离
        parser = PDB.PDBParser(QUIET=True)
        structure = parser.get_structure('complex', pdb_file)
        cys_atoms = {chain1_id: [], chain2_id: []}
        for model in structure:
            for chain in model:
                if chain.id not in [chain1_id, chain2_id]:
                    continue
                for residue in chain:
                    if residue.get_resname() == 'CYS' and 'SG' in residue:
                        cys_atoms[chain.id].append((residue.id[1], residue['SG']))

        print("  所有跨链半胱氨酸对的距离：")
        disulfides = []
        for res1, atom1 in cys_atoms[chain1_id]:
            for res2, atom2 in cys_atoms[chain2_id]:
                distance = atom1 - atom2
                print(f"    {chain1_id}{res1} - {chain2_id}{res2}: {distance:.3f} Å")
                if 1.0 <= distance <= distance_threshold:
                    disulfides.append((res1, res2))

        return {
            "disulfides": disulfides,
            "pdb_file": pdb_file,
            "sequences": {"chain1": seq1, "chain2": seq2}
        }

def main():
    import argparse
    parser = argparse.ArgumentParser(description='预测两条链之间的二硫键（支持批量）')
    parser.add_argument('input_file', help='输入文件：单FASTA文件（含两条序列）或批量格式文件（多块）')
    parser.add_argument('--batch', action='store_true', help='输入为批量格式（标识符+FASTA块）')
    parser.add_argument('--paired', action='store_true', help='输入为多序列FASTA，按顺序两两配对')
    parser.add_argument('--threshold', type=float, default=2.5, help='距离阈值（Å），默认2.5')
    parser.add_argument('--recycle', type=int, default=3, help='迭代次数，默认3')
    parser.add_argument('--outdir', default='results', help='输出目录，默认results')
    args = parser.parse_args()

    os.makedirs(args.outdir, exist_ok=True)

    if args.batch:
        # 批量模式：解析多块文件
        cases = parse_batch_file(args.input_file)
        print(f"共读取到 {len(cases)} 个测试用例")
        for idx, (identifier, seq1, seq2) in enumerate(cases):
            out_json = os.path.join(args.outdir, f"{identifier}.json")
            out_pdb = os.path.join(args.outdir, f"{identifier}.pdb")
            # 断点续跑：如果JSON已存在则跳过
            if os.path.exists(out_json):
                print(f"跳过 {identifier}，结果文件已存在")
                continue
            print(f"\n===== 处理用例 {idx+1}/{len(cases)} : {identifier} =====")
            result = predict_interchain_disulfide(seq1, seq2, num_recycle=args.recycle, distance_threshold=args.threshold, save_pdb_path=out_pdb)
            with open(out_json, 'w') as f:
                json.dump(result, f, indent=2)
            print(f"\n结果已保存至 {out_json}")
            print(f"二硫键: {result['disulfides']}")

    elif args.paired:
        # 成对模式：直接读取FASTA，按顺序两两配对
        records = list(SeqIO.parse(args.input_file, 'fasta'))
        if len(records) % 2 != 0:
            print(f"警告：序列总数为奇数（{len(records)}），最后一个序列将被忽略")
            records = records[:len(records)-1]
        cases = []
        for i in range(0, len(records), 2):
            rec1, rec2 = records[i], records[i+1]
            id1 = rec1.id
            # 提取标识符：尝试匹配 "pair_xxxx"，否则用序号
            match = re.search(r'(pair_\d+)', id1)
            if match:
                pid = match.group(1)
            else:
                pid = f"pair_{i//2+1:04d}"
            cases.append((pid, str(rec1.seq), str(rec2.seq)))
        print(f"共读取到 {len(cases)} 个配对")
        for idx, (identifier, seq1, seq2) in enumerate(cases):
            out_json = os.path.join(args.outdir, f"{identifier}.json")
            out_pdb = os.path.join(args.outdir, f"{identifier}.pdb")
            if os.path.exists(out_json):
                print(f"跳过 {identifier}，结果文件已存在")
                continue
            print(f"\n===== 处理配对 {idx+1}/{len(cases)} : {identifier} =====")
            result = predict_interchain_disulfide(seq1, seq2, num_recycle=args.recycle, distance_threshold=args.threshold, save_pdb_path=out_pdb)
            with open(out_json, 'w') as f:
                json.dump(result, f, indent=2)
            print(f"\n结果已保存至 {out_json}")
            print(f"二硫键: {result['disulfides']}")

    else:
        # 单FASTA模式（两条序列）
        seqs = []
        with open(args.input_file, 'r') as f:
            for record in SeqIO.parse(f, 'fasta'):
                seqs.append(str(record.seq))
        if len(seqs) < 2:
            print("错误：单FASTA模式需要至少两条序列")
            sys.exit(1)
        base = os.path.splitext(os.path.basename(args.input_file))[0]
        out_json = os.path.join(args.outdir, base + '.json')
        out_pdb = os.path.join(args.outdir, base + '.pdb')
        result = predict_interchain_disulfide(seqs[0], seqs[1], num_recycle=args.recycle, distance_threshold=args.threshold, save_pdb_path=out_pdb)
        with open(out_json, 'w') as f:
            json.dump(result, f, indent=2)
        print("\n=== 预测结果 ===")
        print(f"二硫键: {result['disulfides']}")
        print(f"结果已保存至 {out_json}")
        print(f"PDB已保存至 {out_pdb}")

if __name__ == "__main__":
    main()