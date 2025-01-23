#!/bin/bash 
 
 
 
THREADS=1  
REMOVE_TMP=fasle  
FORM="fa" 

while getopts ":i:o:f:t:" opt; do
  case $opt in
    i) INPUT="$OPTARG"
    ;;
    o) OUTPUT_FILE="$OPTARG"
    ;;
    f) FORM="$OPTARG"
    ;;
    t) THREADS="$OPTARG"
    ;; 
    r) REMOVE_TMP=ture 
    ;;
    \?) echo "Invalid option: -$OPTARG" >&2
    ;;
  esac
done
#指定文件位置
SCRIPT_PATH="$(readlink -f "$0")"
BASE_DIR="$(dirname "$SCRIPT_PATH")"
#BASE_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)" 
HMM="$BASE_DIR/hmm/Mn_combine.hmm" 
THREADSHOLD_SCORE="$BASE_DIR/hmm/threshold_score.txt"
 
 
if [ -z "$INPUT" ] || [ -z "$OUTPUT_FILE" ] || [ -z "$THREADS" ] || [ -z "$FORM" ]; then 
  echo "Usage: $0 -i input.fa -o output_file -f {fa|faa} -t threads" 
  exit 1 
fi 

mkdir -p "$OUTPUT_FILE" 
 
# 创建临时文件夹 
TEMP_DIR="./temp_dir_$(date +%s)" 
mkdir "$TEMP_DIR" || exit 
echo "Temporary files will be stored in: $TEMP_DIR" 
TEMP_DIR=$(realpath "$TEMP_DIR")  
#根据-f调整输入 
if [ "$FORM" = "fa" ]
then  
INPUT_BASE=$(basename "$INPUT")  
INPUT_NAME="${INPUT_BASE%.*}" 
prodigal -i $INPUT -a "$TEMP_DIR/${INPUT_NAME}.faa" 
INPUT="$TEMP_DIR/${INPUT_NAME}.faa"  
fi  
if [ "$FORM" = "faa" ] 
then
INPUT_FA="$(realpath "$INPUT")"   
INPUT_BASE=$(basename "$INPUT_FA") 
INPUT_NAME="${INPUT_BASE%.*}" 
fi 




hmmsearch --cpu $THREADS -E 0.01 --tblout "$TEMP_DIR/${INPUT_NAME}_HMM_out" "$BASE_DIR/hmm/Mn_combine.hmm" $INPUT

grep  -v '^#' "$TEMP_DIR/${INPUT_NAME}_HMM_out"  | tr -s ' ' | cut -d" " -f1,3,6 | sed 's/ /\t/g' |sed 's/aln_//' | awk -F'\t' 'NR==FNR{a[$1]=$2;next} {if($3>=a[$2]) print $0}' "$BASE_DIR/hmm/threshold_score.txt" - > "$TEMP_DIR/${INPUT_NAME}_HMM_result"


if [ "$FORM" = "faa" ]
then  
sort -k1,1 "$TEMP_DIR/${INPUT_NAME}_HMM_result" > "$OUTPUT_FILE/${INPUT_NAME}_HMM_result" 
fi 
 
if [ "$FORM" = "fa" ]
then 
 
awk ' 
{ 
    # 使用下划线作为字段分隔符，分割第一列 
    n = split($1, a, "_") 
     
    # 重新组合第一列，去掉最后一个下划线及其后的部分 
    new_first_column = "" 
    for (i = 1; i < n; i++) { 
        new_first_column = new_first_column (i == 1 ? "" : "_") a[i] 
    } 
     
    # 输出修改后的行 
    print new_first_column, $2, $3  # 假设文件有三列，根据实际情况调整 
}' "$TEMP_DIR/${INPUT_NAME}_HMM_result" > "$OUTPUT_FILE/${INPUT_NAME}_HMM_result" 
sort -k 1,1 "$OUTPUT_FILE/${INPUT_NAME}_HMM_result" > ${INPUT_NAME}_tmp && mv ${INPUT_NAME}_tmp "$OUTPUT_FILE/${INPUT_NAME}_HMM_result"
fi 
 
if [ "$REMOVE_TMP" = true ]  
then  
  rm -f $TEMP_DIR  # 删除过程文件的逻辑  
fi 
