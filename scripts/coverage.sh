#!/bin/zsh
CURDIR=`dirname ${BASH_SOURCE[0]-$0}`
LLVM_PROFILE_FILE="$CURDIR/../build/simplicial_arrangement.profraw" $CURDIR/../build/simplicial_arrangement_tests
/usr/local/opt/llvm/bin/llvm-profdata merge -sparse $CURDIR/../build/simplicial_arrangement.profraw -o $CURDIR/../build/simplicial_arrangement.profdata
/usr/local/opt/llvm/bin/llvm-cov report $CURDIR/../build/simplicial_arrangement_tests -instr-profile=$CURDIR/../build/simplicial_arrangement.profdata --ignore-filename-regex="_deps" --ignore-filename-regex="tests/"
/usr/local/opt/llvm/bin/llvm-cov show $CURDIR/../build/simplicial_arrangement_tests -instr-profile=$CURDIR/../build/simplicial_arrangement.profdata --ignore-filename-regex="_deps" --ignore-filename-regex="tests/"  $CURDIR/../src/* -use-color --format html > $CURDIR/../build/coverage.html

