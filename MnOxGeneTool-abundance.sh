#!/bin/bash

# 处理命令行参数
 
THREADS=1 
REMOVE_TMP=fasle 
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
    r) REMOVE_TMP=ture 
    ;;
    \?) echo "Invalid option: -$OPTARG" >&2
    ;;
  esac
done

# 检查必需参数
if [ -z "$INPUT_FA" ] || [ -z "$OUTPUT_FILE" ] || [ -z "$THREADS" ] || [ -z "$ABUNDANCE_MODULE" ]; then
  echo "Usage: $0 -i input.fa -o output_file -a {cell|16s} -t threads"
  exit 1
fi
SCRIPT_PATH="$(readlink -f "$0")"
BASE_DIR="$(dirname "$SCRIPT_PATH")"
#BASE_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
DIAMOND_DB="$BASE_DIR/database/diamonddb/mndb.dmnd"
MNDB_BLAST_DB="$BASE_DIR/database/MOPD/mndb"
SEQUENCE_NAME_FILE="$BASE_DIR/database/MOPD/sequence_name.txt"
LENGTH_FILE="$BASE_DIR/database/MOPD/length.txt"
UNIQUE_16S_FASTA="$BASE_DIR/database/16S/uniques_253bp.fasta"
DB_16S="$BASE_DIR/database/16S/16s_db"
SCG_DB="$BASE_DIR/database/SCG"
SCG_MODEL_DIR="$BASE_DIR/database/SCG_model"

mkdir -p "$OUTPUT_FILE"

# 创建临时文件夹
TEMP_DIR="./temp_dir_$(date +%s)"
mkdir "$TEMP_DIR" || exit
echo "Temporary files will be stored in: $TEMP_DIR"
TEMP_DIR=$(realpath "$TEMP_DIR") 
if [ -n "$SAMPLE" ]; 
then 
INPUT_FA_BASE=$(basename "$INPUT_FA")  
INPUT_FA_NAME="${INPUT_FA_BASE%.*}"
seqkit sample -s 12345 -n "$SAMPLE" "$INPUT_FA" > "$TEMP_DIR/input_sampled.fa" 
INPUT_FA="$TEMP_DIR/input_sampled.fa" 
 
else
INPUT_FA="$(realpath "$INPUT_FA")"  
INPUT_FA_BASE=$(basename "$INPUT_FA") 
INPUT_FA_NAME="${INPUT_FA_BASE%.*}"
fi

# 锰氧化基因计数模块
cd "$TEMP_DIR" || exit



# 执行预筛选
diamond blastx --db "$DIAMOND_DB" --query "$INPUT_FA" --out input_filter_out --outfmt 6 qseqid full_qseq --evalue 10 --id 60 --max-hsps 1 --max-target-seqs 1 --threads "$THREADS" --quiet
# 输出格式转换
cat input_filter_out | while read readline; do
  identifier=$(echo "$readline" | awk -F'\t' '{print $1}')
  sequence=$(echo "$readline" | awk -F'\t' '{$1=""; sub(/^\s+/, ""); print}')
  echo ">$identifier" >> input_filter.fa
  echo "$sequence" >> input_filter.fa
done
# 执行blast
blastx -query input_filter.fa -db "$MNDB_BLAST_DB" -out input_blast_out -evalue 1e-6 -num_threads "$THREADS" -mt_mode 1 -outfmt "6 qseqid qlen sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore"
# 根据identity筛选
awk 'BEGIN {FS="\t"; OFS="\t"} {if (!seen[$1]++ && $4 >= 75) print $0}' input_blast_out > clean_input_out
# 格式整理
awk -F'\t' '{if ($3 ~ /\|/) {split($3, a, "|"); $3 = a[2];} print}' OFS='\t' clean_input_out > clean_input_out_tmp && mv clean_input_out_tmp clean_input_out
# 与数据库建立映射
awk 'NR==FNR{a[$1]=$2;next} $3 in a{print $3, a[$3]}' "$SEQUENCE_NAME_FILE" clean_input_out > mapping.txt
# 拆分每种蛋白的结果
awk 'NR==FNR{map[$1]=$2; next} {output_file = sprintf("output_%s.txt", map[$3]); print > output_file}' mapping.txt clean_input_out
# 创建一个空的文件来存储最终结果
> input_summary.txt
# 遍历所有output_*.txt文件并整合结果


  
# 使用 find 查找文件并处理 
find . -maxdepth 1 -type f -name 'output_*.txt' | while read -r file; do 
  found_files=true 
  annotation=$(basename "$file" | sed -r 's/output_(.+)\.txt/\1/') 
  result=$(sort -k3,3 "$file" | join -t$'\t' -1 3 -2 1 -o 1.1,1.2,1.3,1.4,1.5,2.2 - "$LENGTH_FILE" | awk -F'\t' '{$7=$5/$6; print $7}' OFS='\t' | awk -F'\t' '{sum+=$1} END{print sum}') 
  echo "$annotation $result" >> input_summary.txt
done  
 
if [ ! -s "input_summary.txt" ]; then 
  # 文件为空，输出警告信息 
  echo -e "NA\t0" > input_summary.txt 
  echo "Warning: no mn_oxidizing gene was found" >&2
fi


# 16s计算模块（根据ABUNDANCE_MODULE决定是否执行）
if [ "$ABUNDANCE_MODULE" != "cell" ]; then
  # bwa预过滤
  bwa mem -t "$THREADS" -o input_16s_mem "$UNIQUE_16S_FASTA" "$INPUT_FA"
  # 格式转换
  samtools fasta -F 2308 input_16s_mem > input_16s_fa
  # blast比对
  blastn -db "$DB_16S" -query input_16s_fa -out input_16s_out -evalue 1 -max_hsps 1 -max_target_seqs 1 -mt_mode 1 -outfmt "6 qseqid sseqid pident length slen mismatch gapopen qstart qend sstart send evalue bitscore" -num_threads "$THREADS"
  # 结果格式整理
  awk 'BEGIN {FS="\t"; OFS="\t"} {if (!seen[$1]++) print $0}' input_16s_out > clean_16s_out
  awk '{sum+=$4/253} END {print sum > "input_16s_count"}' clean_16s_out
  if [[ -f "input_summary.txt" && -f "input_16s_count" ]]
then
	N_16s=$(cat "input_16s_count")  
	if [[ -n "$N_16s" && "$N_16s" != 0 ]]; then  
		awk -v N_16s="$N_16s" 'BEGIN {OFS="\t"} {print $1, $2 / N_16s}' "input_summary.txt" > "${OUTPUT_FILE}/${INPUT_FA_NAME}_abundance_16s"  
		else  
		echo "Error: N_16s value is empty or zero"  
	fi  
		else  
		echo "Error: Files input_summary.txt and/or N_16s do not exist"  
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
    if (array_score[$1]) {
      if ($5 > array_score[$1]) {
        array_score[$1] = $5
        array_line[$1, $2] = $4
      }
    } else {
      array_score[$1] = $5
      array_line[$1, $3] = $2  # 修正了原脚本中的错误：$3zz 应为 $3
    }
  } END {
    for (combined in array_line) {
      split(combined, separate, SUBSEP)
      array_length[separate[2]]= array_length[separate[2]] + array_line[combined]
    }
    for (c in array_length) {
      printf "%s\t%s\n", c, array_length[c]
    }
  }' | sort > single_copy_count
  join -1 1 -2 1 single_copy_count "$SCG_DB"/length | awk '{print $1, $2/$3}' > cell_count.txt && rm single_copy_count
  awk '{sum += $2} END {print sum/35}' cell_count.txt > N_cell && rm cell_count.txt
fi

# 根据N_cell和input_summary.txt生成最终输出
if [ -f "input_summary.txt" ] && [ -f "N_cell" ]; then
  N_cell=$(cat "N_cell")
  if [ -n "$N_cell" ] && [ "$N_cell" != 0 ]; then
    awk -v n_cell="$N_cell" 'BEGIN {OFS="\t"} {print $1, $2 / n_cell}' "input_summary.txt" > "${OUTPUT_FILE}/${INPUT_FA_NAME}_abundance_cell"
  else
    echo "Error: N_cell value is empty or zero"
  fi
else
  echo "Error: Files input_summary.txt and/or N_cell do not exist"]
fi 
 
 
#移除过程文件 
if [ "$REMOVE_TMP" = true ] 
then 
  rm -f $TEMP_DIR  # 删除过程文件的逻辑 
fi 

