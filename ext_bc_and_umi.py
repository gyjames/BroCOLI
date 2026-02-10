# 2026.01.23
# panwei
# 从bam文件中提取barcode和umi, 这里是cellreanger的, 标签是CB:Z:和UB:Z:
# python ext_bc_and_umi.py --bam possorted_genome_bam.bam -f filter_barcodes_clean.tsv -o sc_bc_umi.txt

import argparse
import pysam

# 创建参数解析器
parser = argparse.ArgumentParser(description="generate reference sequence!")
# 添加参数
parser.add_argument("-b", "--bam", required=True, help="input bam file")
parser.add_argument("-f", "--filter", default=None, help="filter barcode tsv file")
parser.add_argument("-o", "--output", required=True, help="output txt path")

# 解析参数
args = parser.parse_args()

# 1. 设置输入和输出文件路径
bam_file_path = args.bam
filter_barcode_path = args.filter
output_temp_file = args.output

white_set = set()
if len(filter_barcode_path) != 0:
    with open(args.filter, "r") as infile:
        for line in infile:
            cols = line.rstrip("\n")
            white_set.add(cols)


def extract_barcode_umi():
    bam = pysam.AlignmentFile(bam_file_path, "rb")
    
    with open(output_temp_file, "w") as f_out:
        print(f"Processing {bam_file_path}, please wait...")
        
        count = 0
        # 逐行读取 reads
        for read in bam:
            # 只有当这条 read 既有 CB 又有 UB 标签时才提取
            if read.has_tag("CB") and read.has_tag("UB"):
                # get_tag 自动只取值，不带 CB:Z: 前缀
                # cb = read.get_tag("CB")
                cb = read.get_tag("CB").split('-')[0]
                ub = read.get_tag("UB")
                
                if len(white_set) > 0:
                    if cb in white_set:
                        f_out.write(f"{cb}\t{ub}\n")
                else:
                    f_out.write(f"{cb}\t{ub}\n")
                
            count += 1
            if count % 1000000 == 0:
                print(f"Processing {count} reads...")
                
    bam.close()
    print("Extract completed!")

if __name__ == "__main__":
    extract_barcode_umi()


