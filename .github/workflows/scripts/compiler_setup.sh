#/bin/bash

echo $CC 
echo $CXX

export TOOLCHAIN=$1

sudo apt install ${CC} ${CXX} ${FC}

sed -i -e "s/CI_CXX/${CXX}/g"  ${TOOLCHAIN}
sed -i -e "s/CI_C/${CC}/g" ${TOOLCHAIN}
sed -i -e "s/CI_FC/${FC}/g" ${TOOLCHAIN}
echo $PWD && cat ${TOOLCHAIN}
