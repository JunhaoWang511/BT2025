#!/bin/bash
# 这个脚本是对UI服务器上/ustcfs3/stcf/BT2025/目录下12月的束流数据进行处理
para_num() {
    if test $# -lt 1; then
        echo "请输入至少2个参数, par1: UI服务器上数据文件夹的绝对路径; par2: 处理文件个数"
        exit 1
    fi
}
get_last_dirs() {
    local path="$1"
    local n="$2"
    IFS='/' read -r -a parts <<<"$path"
    local total=${#parts[@]}
    for ((i = total - n; i < total; i++)); do
        eval "dir$((i - total + n + 1))=${parts[$i]}"
    done
}
# ！！输入参数1: UI01服务器中的数据文件夹绝对路径“**/ECAL/”
# ！！输入参数2：处理的文件数目(缺省时拷贝全部文件)
DATAFILEPATH=$1
FILENUMBER=$2
para_num $@

DATAFILEPATH=${DATAFILEPATH%/}
if [ $(basename $DATAFILEPATH) != "ECAL" ]; then
    DATAFILEPATH=$DATAFILEPATH/ECAL
fi
# 获取最后四级目录名到dir1、2、3、4
get_last_dirs $DATAFILEPATH 4
# 脚本所在路径
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
DATA_DIR=${SCRIPT_DIR}/data/"$dir1$dir2$dir3"
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
    i=0
    # 解码前FILENUMBER个文件
    realpath $DATAFILEPATH/data_ECAL*.dat >list.txt
    while read -r temp; do
        dataname=$(basename $temp)
        if [[ "$dataname" == data*ECAL*.dat ]]; then
            # 后台解码文件
            ${DECODE_DIR}/ECALdig2root2025 ${temp} ${dataname%.dat}_decode.root >/dev/null 2>&1 &
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
# 重建在线提取数据
if [ ! -f "rec_online.root" ]; then
    realpath data*ECAL*decode.root >list.txt
    while read -r temp; do
        dataname=$(basename $temp)
        ${DECODE_DIR}/multiReconstruction ${dataname} ${dataname%_decode.root}_rec_online.root >/dev/null 2>&1 &
    done <list.txt
    rm list.txt
    # 等待后台程序执行完毕
    echo "wait onlineTQ reconstruction processing..."
    wait
    hadd rec_online.root data*ECAL*rec_online.root
    echo "generate $(pwd)/rec_online.root"
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
echo "运行结束，耗时 $((SECONDS/60)) 分 $((SECONDS%60)) 秒"
