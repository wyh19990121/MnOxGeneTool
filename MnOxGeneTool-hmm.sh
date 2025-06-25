#!/bin/bash 

# 使用环境变量定位资源文件
if [ -z "$MNOX_DATA_DIR" ]; then
    echo "ERROR: MNOX_DATA_DIR environment variable not set."
    echo "Please activate the conda environment for this tool."
    exit 1
fi

DATA_DIR="$MNOX_DATA_DIR"  # 使用环境变量定义的数据目录

THREADS=1  
REMOVE_TMP=false  
FORM="fa" 

# 添加 'r' 选项到 getopts
while getopts ":i:o:f:t:r" opt; do
  case $opt in
    i) INPUT="$OPTARG"
    ;;
    o) OUTPUT_FILE="$OPTARG"
    ;;
    f) FORM="$OPTARG"
    ;;
    t) THREADS="$OPTARG"
    ;; 
    r) REMOVE_TMP=true  
    ;;
    \?) echo "Invalid option: -$OPTARG" >&2
    exit 1
    ;;
  esac
done

# 使用 DATA_DIR 定义资源路径
HMM="$DATA_DIR/hmm/Mn_combine.hmm" 
THREADSHOLD_SCORE="$DATA_DIR/hmm/threshold_score.txt"

# 验证资源文件存在
if [ ! -f "$HMM" ]; then
    echo "ERROR: HMM file not found at $HMM"
    exit 1
fi

if [ ! -f "$THREADSHOLD_SCORE" ]; then
    echo "ERROR: Threshold score file not found at $THREADSHOLD_SCORE"
    exit 1
fi

# 输入验证
if [ -z "$INPUT" ] || [ -z "$OUTPUT_FILE" ] || [ -z "$THREADS" ] || [ -z "$FORM" ]; then 
  echo "Usage: $0 -i input.fa -o output_file -f {fa|faa} -t threads [-r]" 
  exit 1 
fi 

mkdir -p "$OUTPUT_FILE" 

# 创建临时文件夹 
TEMP_DIR="$OUTPUT_FILE/temp_dir_$(date +%s)"  # 将临时目录放在输出目录下
mkdir "$TEMP_DIR" || exit 
echo "Temporary files will be stored in: $TEMP_DIR" 

# 根据-f调整输入 
if [ "$FORM" = "fa" ]
then  
    INPUT_BASE=$(basename "$INPUT")  
    INPUT_NAME="${INPUT_BASE%.*}" 
    prodigal -i $INPUT -a "$TEMP_DIR/${INPUT_NAME}.faa" 
    INPUT="$TEMP_DIR/${INPUT_NAME}.faa"  
elif [ "$FORM" = "faa" ] 
then
    INPUT_FA="$(realpath "$INPUT")"   
    INPUT_BASE=$(basename "$INPUT_FA") 
    INPUT_NAME="${INPUT_BASE%.*}" 
else
    echo "Invalid format: $FORM. Use 'fa' or 'faa'."
    exit 1
fi 

# 使用绝对路径运行 hmmsearch
hmmsearch --cpu $THREADS -E 0.01 --tblout "$TEMP_DIR/${INPUT_NAME}_HMM_out" "$HMM" "$INPUT"

# 处理结果
grep -v '^#' "$TEMP_DIR/${INPUT_NAME}_HMM_out" | tr -s ' ' | cut -d" " -f1,3,6 | sed 's/ /\t/g' | sed 's/aln_//' | awk -F'\t' 'NR==FNR{a[$1]=$2;next} {if($3>=a[$2]) print $0}' "$THREADSHOLD_SCORE" - > "$TEMP_DIR/${INPUT_NAME}_HMM_result"

# 根据格式处理结果
if [ "$FORM" = "faa" ]
then  
    sort -k1,1 "$TEMP_DIR/${INPUT_NAME}_HMM_result" > "$OUTPUT_FILE/${INPUT_NAME}_HMM_result" 
elif [ "$FORM" = "fa" ]
then 
    awk '{
        n = split($1, a, "_")
        new_first_column = ""
        for (i = 1; i < n; i++) {
            new_first_column = new_first_column (i == 1 ? "" : "_") a[i]
        }
        print new_first_column, $2, $3
    }' "$TEMP_DIR/${INPUT_NAME}_HMM_result" > "$OUTPUT_FILE/${INPUT_NAME}_HMM_result" 
    sort -k1,1 "$OUTPUT_FILE/${INPUT_NAME}_HMM_result" -o "$OUTPUT_FILE/${INPUT_NAME}_HMM_result"
fi 

# 清理临时文件
if [ "$REMOVE_TMP" = true ]  
then  
    echo "Cleaning temporary files..."
    rm -rf "$TEMP_DIR"
else
    echo "Temporary files kept at: $TEMP_DIR"
fi

echo "Analysis completed. Results saved to: $OUTPUT_FILE"
