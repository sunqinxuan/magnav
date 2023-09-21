# magnav

## requirements

HDF5 package is required: [HDF5](https://portal.hdfgroup.org/display/HDF5/HDF5)

tested on Ubuntu 20.04.

## build 

```sh
git clone https://github.com/sunqinxuan/magnav.git
cd magnav
mkdir build
cd build 
cmake ..
make 
```

## run

### dataset

Download the dataset from [here](https://magnav.mit.edu/) and put them into a `data` directory.

Remember to change the dataset path in `src/magnav.cpp` file.

### run TL-model 

```sh
bin/magnav
```

A file `data_TL.h5` is generated by this process.

### run ML-model

Run `untitled.m` in Matlab.

