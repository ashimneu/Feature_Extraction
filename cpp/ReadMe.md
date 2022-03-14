To compile and build executable:

1. run prepare_cmake.sh using bash command:
```
$./prepare_cmake.sh [build_directory]
```

2. This will create a subdirectory with the name build_directory. If no name is provided, the default name is build

3. change to the build directory and make the project
```
$cd build
$make -jN
```
N should be substited with the number of threads we want to use to make the project. The more threads dedicated, the faster the build, defaults to 4.


4. run program. When cmake has finished building the project, there will be an executable named feature_detect
```
$./feature_detect -u [input_pcd_file_name]
```
This will load the pcd file listed. The file should be located in the same directory as the executable.
-u command line argument is used to indicate the subsequent string is the pcd file name
