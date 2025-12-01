#!/bin/bash
# 这个脚本是对UI服务器上/ustcfs3/stcf/BT2025/目录下的束流数据进行处理
para_num() {
    if test $# -lt 2; then
        echo "请输入至少2个参数, par1: DAQ服务器上数据文件夹的绝对路径;par2: 处理文件个数"
        exit 1
    fi
}
# ！！输入参数1: UI01服务器中的数据文件夹绝对路径“**/ECAL/”
# ！！输入参数2：处理的文件数目(缺省时拷贝全部文件)
DATAFILEPATH=$1
FILENUMBER=$2
para_num $@

DATERUN=$(basename $DATAFILEPATH)
# 脚本所在路径
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
DATA_DIR=${SCRIPT_DIR}/data/"$DATERUN"
mkdir -p $DATA_DIR
echo "make data directory $(realpath $DATA_DIR)"
pushd ${DATA_DIR} >/dev/null

# 数据处理程序路径
DECODE_DIR=${SCRIPT_DIR}/../../build

# 解码文件生成decode.root
if [ ! -f "decode.root" ]; then
    # 判断是否生成解码程序
    if [ ! -f "${DECODE_DIR}/ECALdig2root2025" ]; then
        echo "Decode programe dose not exist!"
        exit 1
    fi
    # 解码前FILENUMBER个文件
    if [ -n "$FILENUMBER" ]; then
        ls $DATAFILEPATH/data_ECAL*.dat | sort | head -n "$FILENUMBER" | while read -r f; do
            datapath=$(realpath "$f")
            ${DECODE_DIR}/ECALdig2root2025 ${datapath} "$(basename "${datapath%.dat}").root"
        done
    else
        ls $DATAFILEPATH/data_ECAL*.dat | sort | while read -r f; do
            datapath=$(realpath "$f")
            ${DECODE_DIR}/ECALdig2root2025 ${datapath} "$(basename "${datapath%.dat}").root"
        done
    fi
    # 等待后台程序执行完毕
    # wait
    hadd decode.root data*ECAL*.root
    echo "generate $(pwd)/decode.root"
    # 移除中间数据文件
    rm data*ECAL*
fi

# 数字化和重建
if [ ! -f 'digi.root' ]; then
    ${DECODE_DIR}/ECALDigi $(pwd)/decode.root digi.root 
    echo "generate $(pwd)/digi.root"
fi
if [ ! -f "rec.root" ]; then
    ${DECODE_DIR}/Reconstruction $(pwd)/digi.root rec.root
    echo "generate $(pwd)/rec.root"
fi

# echo "fit electron energy spectrum"
# ${SCRIPT_DIR}/DrawEnergy rec.root

# echo "fit electron position distribution"
# ${SCRIPT_DIR}/DrawPosition rec.root Tracker-step4-rec.root

# echo "print working state information"
# ${SCRIPT_DIR}/Draw5x5NoiseTempGratio decode.root

popd >/dev/null
