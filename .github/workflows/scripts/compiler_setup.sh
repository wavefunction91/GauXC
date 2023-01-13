#/bin/bash

echo $CC 
echo $CXX

export TOOLCHAIN=$1

if [[ "$CC" == "clang"* ]]
then
  sudo apt install ${CC}
else
  sudo apt install ${CC} ${CXX}
fi

sed -i -e "s/CI_CXX/${CXX}/g"  ${TOOLCHAIN}
sed -i -e "s/CI_C/${CC}/g" ${TOOLCHAIN}
if [[ "$CC" == "clang"* ]]
then
  echo "set( GAUXC_ENABLE_OPENMP OFF CACHE BOOL \"\" FORCE )" >> ${TOOLCHAIN}
fi
echo $PWD && cat ${TOOLCHAIN}
