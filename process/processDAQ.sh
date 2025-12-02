#!/bin/bash
# 这个脚本是对DAQ服务器上lunon下的束流数据进行处理
para_num() {
  if test $# -lt 1; then
    echo "Please Input data file realpath(/data/DAQ/Beam_2510/*.dat), or plus 2nd parameter: file number to copy"
    exit 1
  fi
}
para_num $@

# 输入参数1: DAQ服务器中的数据文件夹路径or'202510xx_xx'
# 输入参数2：拷贝的文件数目(缺省时拷贝全部文件)
# 输入参数3：电子束流能量(缺省时自动寻找能量)
DATAFILEPATH=$1
FILENUMBER=$2
ENERGYCUT=$3

# 数据文件夹（日期）名
DATENAME=$(basename $DATAFILEPATH)

# 脚本所在路径
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

mkdir -p ${SCRIPT_DIR}/data/${DATENAME}
pushd ${SCRIPT_DIR}/data/${DATENAME} >/dev/null

# 数据处理程序路径
DECODE_DIR=${SCRIPT_DIR}/../../build

# 解码文件生成decode.root
if [ ! -f "decode.root" ]; then
  # 判断是否生成解码程序
  if [ ! -f "${DECODE_DIR}/ECALdig2root2025" ]; then
    echo "Decode programe dose not exist!"
    exit 1
  fi
  i=0
  if [[ ! "$DATAFILEPATH" == *"/"* ]]; then
    DATAFILEPATH=/data/DAQ/Beam_2510/$DATAFILEPATH
  fi
  # 拷贝DAQ服务器的数据到本地，跳过已有文件
  ssh "lunon@192.168.7.150" "find '$DATAFILEPATH' -type f" >list.txt
  while read -r temp; do
    if [[ "$(basename $temp)" == data*ECAL*.dat ]]; then
      # 后台拷贝文件和解码
      dataname=$(basename $temp)
      rsync -h --ignore-existing --progress lunon@192.168.7.150:${DATAFILEPATH}/${dataname} .
      ${DECODE_DIR}/ECALdig2root2025 ${dataname} ${dataname%.dat}_decode.root >/dev/null 2>&1 &
      ((i++))
      if [[ (-n $FILENUMBER) && ($i -eq $FILENUMBER) ]]; then
        break
      fi
    fi
  done <list.txt
  rm list.txt
  # 等待后台程序执行完毕
  echo "wait decode processing..."
  wait
  hadd decode.root data*ECAL*.root
  echo "generate $(pwd)/decode.root"
fi

# 数字化和重建
if [ ! -f 'digi.root' ]; then
  realpath data*ECAL*decode.root >list.txt
  while read -r temp; do
    dataname=$(basename $temp)
    ${DECODE_DIR}/ECALDigi ${dataname} ${dataname%_decode.root}_digi.root >/dev/null 2>&1 &
  done <list.txt
  rm list.txt
  # 等待后台程序执行完毕
  echo "wait digitization processing..."
  wait
  hadd digi.root data*ECAL*digi.root
  echo "generate $(pwd)/digi.root"
fi
if [ ! -f "rec.root" ]; then
  realpath data*ECAL*digi.root >list.txt
  while read -r temp; do
    dataname=$(basename $temp)
    ${DECODE_DIR}/Reconstruction ${dataname} ${dataname%_digi.root}_rec.root >/dev/null 2>&1 &
  done <list.txt
  rm list.txt
  # 等待后台程序执行完毕
  echo "wait reconstruction processing..."
  wait
  hadd rec.root data*ECAL*rec.root
  echo "generate $(pwd)/rec.root"
fi
# 移除中间数据文件
rm data*ECAL*

# echo "fit electron energy spectrum"
# ${SCRIPT_DIR}/DrawEnergy rec.root

# echo "fit electron position distribution"
# ${SCRIPT_DIR}/DrawPosition rec.root Tracker-step4-rec.root

# echo "print working state information"
# ${SCRIPT_DIR}/Draw5x5NoiseTempGratio decode.root

popd >/dev/null
