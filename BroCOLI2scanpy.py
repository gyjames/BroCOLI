import argparse
import os
import pandas as pd
import scipy.sparse as sp
from scipy.io import mmwrite
import gzip, shutil

# 创建参数解析器
parser = argparse.ArgumentParser(description='Convert BroCOLI matrix to seurat and scanpy!')
# 添加参数
parser.add_argument("-i", "--input", required=True, help="input gene/isoform barcode expression matrix")
parser.add_argument("-g", "--gtf", required=False, default=None, help="gene/isoform name")
parser.add_argument("-t", "--type", required=False, default="gene", help="gene/isoform name")
parser.add_argument("-o", "--output", required=True, help="output folder path")
# 解析参数
args = parser.parse_args()

df = pd.read_csv(args.input, sep='\t', index_col=0)
df.columns = [col[5:] + "-1" for col in df.columns]

if args.type != "gene":
    df = df.drop(columns=["gene_id"])

print("Before the merge:", df.shape)
df = df.groupby(df.index).sum()
print("After the merge:", df.shape)

barcode_output_path = os.path.join(args.output, "barcodes.tsv.gz")
# 1. 保存 barcodes.tsv.gz
df.columns.to_series().to_csv(barcode_output_path, index=False, header=False)

# 2. 保存 features.tsv.gz
# 10x 要求三列：gene_id  gene_name  feature_type
# 如果你只有 transcript_id，就把 gene_id 和 gene_name 都写成 transcript_id
# 如果有gtf, 从gtf中解析已有的转录本及其对应的名称;
if args.type != "gene":
    transcripts = []
    if (len(args.gtf) != 0):
        with open(args.gtf) as f:
            for line in f:
                if line.startswith("#"):
                    continue
                fields = line.strip().split("\t")
                if fields[2] != "transcript":  
                    continue
                attr_field = fields[8]
                attrs = {}
                for attr in attr_field.split(";"):
                    if attr.strip():
                        key, value = attr.strip().split(" ", 1)
                        attrs[key] = value.strip('"')

                transcript_id = attrs.get("transcript_id")
                transcript_name = attrs.get("transcript_name")

                if transcript_id and transcript_name:
                    transcripts.append([transcript_id, transcript_name])
    gtf_df = pd.DataFrame(transcripts, columns=["transcript_id", "transcript_name"])
    features = (
        pd.DataFrame({"gene_id": df.index})
        .merge(gtf_df, left_on="gene_id", right_on="transcript_id", how="left")
        .drop(columns=["transcript_id"])
    )
    features["transcript_name"] = features["transcript_name"].fillna(features["gene_id"])
    features["feature_type"] = "Gene Expression"
    features = features.rename(columns={"transcript_name": "gene_name"})

else:
    genes = []
    if (len(args.gtf) != 0):
        with open(args.gtf) as f:
            for line in f:
                if line.startswith("#"):
                    continue
                fields = line.strip().split("\t")
                if fields[2] != "gene":  
                    continue         
                attr_field = fields[8]
                attrs = {}
                for attr in attr_field.split(";"):
                    if attr.strip():
                        key, value = attr.strip().split(" ", 1)
                        attrs[key] = value.strip('"')
                gene_id = attrs.get("gene_id")
                gene_name = attrs.get("gene_name")
                if gene_id and gene_name:
                    genes.append([gene_id, gene_name])
    gtf_df = pd.DataFrame(genes, columns=["gene_id", "gene_name"])
    features = (
        pd.DataFrame({"gene_id": df.index})
        .merge(gtf_df, on="gene_id", how="left")
    )
    features.insert(2, "source", "Gene Expression")

feature_output_path = os.path.join(args.output, "features.tsv.gz")
features.to_csv(feature_output_path, sep="\t", index=False, header=False)

# 3. 保存 matrix.mtx.gz
# 转换成稀疏矩阵(注意 MatrixMarket 是 (genes × barcodes))
matrix_output_path = os.path.join(args.output, "matrix.mtx.gz")
mat = sp.csr_matrix(df.values)
mmwrite("matrix.mtx", mat)   # 先写成 matrix.mtx
# 压缩
with open("matrix.mtx", "rb") as f_in, gzip.open(matrix_output_path, "wb") as f_out:
    shutil.copyfileobj(f_in, f_out)
# 删除临时文件
os.remove("matrix.mtx")








