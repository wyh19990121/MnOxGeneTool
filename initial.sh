#!/bin/bash



INSTALL_DIR=/usr/local/bin

# 处理命令行选项
while getopts ":i:" opt; do
  case $opt in
    i)
      INSTALL_DIR="$OPTARG"
      ;;
    \?)
      echo "Invalid option: -$OPTARG" >&2
      exit 1
      ;;
  esac
done

# 检查源文件是否存在
if [[ ! -f "MnOxGeneTool-abundance.sh" ]]; then
  echo "Error: MnOxGeneTool-abundance.sh not found!" >&2
  exit 1
fi

if [[ ! -f "MnOxGeneTool-hmm.sh" ]]; then
  echo "Error: MnOxGeneTool-hmm.sh not found!" >&2
  exit 1
fi

# 设置源文件可执行权限
chmod +x MnOxGeneTool-abundance.sh
chmod +x MnOxGeneTool-hmm.sh

# 创建符号链接
ln -sf "$(pwd)/MnOxGeneTool-abundance.sh" "${INSTALL_DIR}/MnOxGeneTool-abundance"
ln -sf "$(pwd)/MnOxGeneTool-hmm.sh" "${INSTALL_DIR}/MnOxGeneTool-hmm"

# 确保链接后的文件也是可执行的
chmod +x "${INSTALL_DIR}/MnOxGeneTool-abundance"
chmod +x "${INSTALL_DIR}/MnOxGeneTool-hmm"

echo "Installation completed successfully. Scripts are installed in ${INSTALL_DIR}."