#!/bin/bash

# 使用环境变量定位资源文件
if [ -z "$MNOX_DATA_DIR" ]; then
    echo "ERROR: MNOX_DATA_DIR environment variable not set."
    echo "Please activate the conda environment for this tool."
    exit 1
fi

DATA_DIR="$MNOX_DATA_DIR"  # 使用环境变量定义的数据目录

# 处理命令行参数
THREADS=1 
REMOVE_TMP=false 
ABUNDANCE_MODULE=1

while getopts ":i:o:a:t:s:r" opt; do
  case $opt in
    i) INPUT_FA="$OPTARG"
    ;;
    o) OUTPUT_FILE="$OPTARG"
    ;;
    a) ABUNDANCE_MODULE="$OPTARG"
    ;;
    t) THREADS="$OPTARG" 
    ;; 
    s) SAMPLE="$OPTARG"
    ;; 
    r) REMOVE_TMP=true  
    ;;
    \?) echo "Invalid option: -$OPTARG" >&2
    exit 1
    ;;
  esac
done

# 检查必需参数
if [ -z "$INPUT_FA" ] || [ -z "$OUTPUT_FILE" ] || [ -z "$THREADS" ] || [ -z "$ABUNDANCE_MODULE" ]; then
  echo "Usage: $0 -i input.fa -o output_file -a {cell|16s} -t threads [-s sample_size] [-r]"
  exit 1
fi

# 使用 DATA_DIR 定义资源路径
DIAMOND_DB="$DATA_DIR/database/diamonddb/mndb.dmnd"
MNDB_BLAST_DB="$DATA_DIR/database/MOPD/mndb"
SEQUENCE_NAME_FILE="$DATA_DIR/database/MOPD/sequence_name.txt"
LENGTH_FILE="$DATA_DIR/database/MOPD/length.txt"
UNIQUE_16S_FASTA="$DATA_DIR/database/16S/uniques_253bp.fasta"
DB_16S="$DATA_DIR/database/16S/16s_db"
SCG_DB="$DATA_DIR/database/SCG"
SCG_MODEL_DIR="$DATA_DIR/database/SCG_model"

# 验证资源文件存在
if [ ! -d "$SCG_DB" ]; then
    echo "ERROR: SCG database not found at $SCG_DB"
    exit 1
fi

if [ ! -d "$SCG_MODEL_DIR" ]; then
    echo "ERROR: SCG model directory not found at $SCG_MODEL_DIR"
    exit 1
fi

mkdir -p "$OUTPUT_FILE"

# 创建临时文件夹（放在输出目录下）
TEMP_DIR="$OUTPUT_FILE/temp_dir_$(date +%s)"
mkdir "$TEMP_DIR" || exit
echo "Temporary files will be stored in: $TEMP_DIR"

# 采样处理
if [ -n "$SAMPLE" ]; then 
    INPUT_FA_BASE=$(basename "$INPUT_FA")  
    INPUT_FA_NAME="${INPUT_FA_BASE%.*}"
    seqkit sample -s 12345 -n "$SAMPLE" "$INPUT_FA" > "$TEMP_DIR/input_sampled.fa" 
    INPUT_FA="$TEMP_DIR/input_sampled.fa" 
else
    INPUT_FA_BASE=$(basename "$INPUT_FA") 
    INPUT_FA_NAME="${INPUT_FA_BASE%.*}"
fi

# 锰氧化基因计数模块
cd "$TEMP_DIR" || exit

# 执行预筛选
diamond blastx --db "$DIAMOND_DB" --query "$INPUT_FA" --out input_filter_out \
    --outfmt 6 qseqid full_qseq --evalue 10 --id 60 --max-hsps 1 \
    --max-target-seqs 1 --threads "$THREADS" --quiet

# 输出格式转换
cat input_filter_out | while read -r readline; do
    identifier=$(echo "$readline" | awk -F'\t' '{print $1}')
    sequence=$(echo "$readline" | awk -F'\t' '{$1=""; sub(/^\s+/, ""); print}')
    echo ">$identifier" >> input_filter.fa
    echo "$sequence" >> input_filter.fa
done

# 执行blast
blastx -query input_filter.fa -db "$MNDB_BLAST_DB" -out input_blast_out \
    -evalue 1e-6 -num_threads "$THREADS" -mt_mode 1 \
    -outfmt "6 qseqid qlen sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore"

# 根据identity筛选
awk 'BEGIN {FS="\t"; OFS="\t"} !seen[$1]++ && $4 >= 75 {print}' input_blast_out > clean_input_out

# 格式整理
awk -F'\t' '{if ($3 ~ /\|/) {split($3, a, "|"); $3 = a[2];} print}' OFS='\t' clean_input_out > clean_input_out.tmp
mv clean_input_out.tmp clean_input_out

# 与数据库建立映射
awk 'NR==FNR{a[$1]=$2;next} $3 in a {print $3, a[$3]}' "$SEQUENCE_NAME_FILE" clean_input_out > mapping.txt

# 拆分每种蛋白的结果
awk 'NR==FNR{map[$1]=$2; next} {
    if ($3 in map) {
        output_file = "output_" map[$3] ".txt"
        print > output_file
    }
}' mapping.txt clean_input_out

# 创建一个空的文件来存储最终结果
> input_summary.txt

# 遍历所有output_*.txt文件并整合结果
found_files=false
for file in output_*.txt; do
    [ -f "$file" ] || continue
    found_files=true
    annotation=$(basename "$file" | sed -E 's/output_(.+)\.txt/\1/')
    result=$(sort -k3,3 "$file" | join -t$'\t' -1 3 -2 1 -o 1.1,1.2,1.3,1.4,1.5,2.2 - "$LENGTH_FILE" | \
        awk -F'\t' '{cov = $5/$6; print cov}' | \
        awk '{sum += $1} END {print sum}')
    echo -e "$annotation\t$result" >> input_summary.txt
done

if [ "$found_files" = false ]; then 
    echo -e "NA\t0" > input_summary.txt 
    echo "Warning: no mn_oxidizing gene was found" >&2
fi

# 16s计算模块（根据ABUNDANCE_MODULE决定是否执行）
if [ "$ABUNDANCE_MODULE" != "cell" ]; then
    # bwa预过滤
    bwa mem -t "$THREADS" -o input_16s_mem.sam "$UNIQUE_16S_FASTA" "$INPUT_FA"
    
    # 格式转换
    samtools fasta -F 2308 input_16s_mem.sam > input_16s_fa.fasta
    
    # blast比对
    blastn -db "$DB_16S" -query input_16s_fa.fasta -out input_16s_out \
        -evalue 1 -max_hsps 1 -max_target_seqs 1 -mt_mode 1 \
        -outfmt "6 qseqid sseqid pident length slen mismatch gapopen qstart qend sstart send evalue bitscore" \
        -num_threads "$THREADS"
    
    # 结果格式整理
    awk 'BEGIN {FS="\t"; OFS="\t"} !seen[$1]++' input_16s_out > clean_16s_out
    
    # 计算16s数量
    awk '{sum += $4/253} END {print sum > "input_16s_count"}' clean_16s_out
    
    # 生成丰度结果
    if [[ -f "input_summary.txt" && -f "input_16s_count" ]]; then
        N_16s=$(cat "input_16s_count")  
        if [[ -n "$N_16s" && "$N_16s" != "0" && "$N_16s" != "0.0" ]]; then  
            awk -v N_16s="$N_16s" 'BEGIN {OFS="\t"} {print $1, $2 / N_16s}' "input_summary.txt" > "${OUTPUT_FILE}/${INPUT_FA_NAME}_abundance_16s"  
        else  
            echo "Warning: N_16s value is zero or invalid" > "${OUTPUT_FILE}/${INPUT_FA_NAME}_abundance_16s"
        fi  
    else  
        echo "Error: Missing input files for 16s calculation" > "${OUTPUT_FILE}/${INPUT_FA_NAME}_abundance_16s"
    fi  
fi

# cellcount计算模块（根据ABUNDANCE_MODULE决定是否执行）
if [ "$ABUNDANCE_MODULE" != "16s" ]; then
    # 预测蛋白
    FragGeneScanRs -s "$INPUT_FA" -a input.faa -w 0 -t illumina_5 -p "$THREADS"
    
    # 蛋白比对到35个单拷贝基因数据库
    uproc-prot --threads "$THREADS" --output input_uprot_out --preds --pthresh 3 "$SCG_DB" "$SCG_MODEL_DIR" input.faa
    
    # 计算细胞数
    cut -f1,3,4,5 -d"," "input_uprot_out" | \
    awk 'BEGIN {OFS="\t"; FS=","} {
        key = $1
        score = $4
        gene = $3
        len = $2
        
        # 如果这个contig还没有记录，或者当前得分更高
        if (!(key in best_score) || (score > best_score[key]) {
            best_score[key] = score
            best_gene[key] = gene
            best_length[key] = len
        }
    } END {
        # 累加每个基因的长度
        for (key in best_gene) {
            gene = best_gene[key]
            gene_lengths[gene] += best_length[key]
        }
        
        # 输出每个基因的总长度
        for (gene in gene_lengths) {
            print gene, gene_lengths[gene]
        }
    }' | sort > single_copy_count
    
    # 计算每个基因的覆盖度
    join -1 1 -2 1 single_copy_count "$SCG_DB/length" | awk '{print $1, $2/$3}' > cell_count.txt
    rm -f single_copy_count
    
    # 计算平均细胞数
    awk '{sum += $2} END {if (NR > 0) print sum/35; else print 0}' cell_count.txt > N_cell
    rm -f cell_count.txt
fi

# 根据N_cell和input_summary.txt生成最终输出
if [ -f "input_summary.txt" ] && [ -f "N_cell" ]; then
    N_cell=$(cat "N_cell")
    if [[ -n "$N_cell" && "$N_cell" != "0" && "$N_cell" != "0.0" ]]; then
        awk -v n_cell="$N_cell" 'BEGIN {OFS="\t"} {print $1, $2 / n_cell}' "input_summary.txt" > "${OUTPUT_FILE}/${INPUT_FA_NAME}_abundance_cell"
    else
        echo "Warning: N_cell value is zero or invalid" > "${OUTPUT_FILE}/${INPUT_FA_NAME}_abundance_cell"
    fi
else
    echo "Error: Missing input files for cell calculation" > "${OUTPUT_FILE}/${INPUT_FA_NAME}_abundance_cell"
fi 
 
# 移除过程文件 
if [ "$REMOVE_TMP" = true ]; then 
    echo "Cleaning temporary files..."
    rm -rf "$TEMP_DIR"
else
    echo "Temporary files kept at: $TEMP_DIR"
fi

echo "Analysis completed. Results saved to: $OUTPUT_FILE"
