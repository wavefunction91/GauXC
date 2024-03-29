#/bin/bash

export CSUITE=$1
export CVER=$2

if [[ "${CSUITE}" == "llvm" ]]
then
  update-alternatives --set clang   /usr/bin/clang-${CVER}
  update-alternatives --set clang++ /usr/bin/clang++-${CVER}
  update-alternatives --install /usr/bin/cc  cc  /usr/bin/clang   30
  update-alternatives --install /usr/bin/c++ c++ /usr/bin/clang++ 30
elif [[ "${CSUITE}" == "gnu" ]]
then
  update-alternatives --set gcc /usr/bin/gcc-${CVER}
  update-alternatives --set g++ /usr/bin/g++-${CVER}
  update-alternatives --install /usr/bin/cc  cc  /usr/bin/gcc 30
  update-alternatives --install /usr/bin/c++ c++ /usr/bin/g++ 30
else
  echo "Compiler Suite Not Recognized!"
  exit 125
fi

