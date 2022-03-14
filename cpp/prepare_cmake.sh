#!/bin/bash
if (( $# < 1 ))
then
    rootdir="build"
else
    rootdir=$1
fi

mkdir -p $rootdir

# sources=( "data" "conf" )

# for dir1 in "${sources[@]}"
# do
#     ln -s ../$dir1 $rootdir/$dir1
#     echo $dir1
# done

cd $rootdir
cmake .. -DCMAKE_EXPORT_COMPILE_COMMANDS=true -Wno-dev
