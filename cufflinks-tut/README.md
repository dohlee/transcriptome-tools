# Cufflinks Tutorial

More extensive documentation can be found at [the official github repository](https://github.com/cole-trapnell-lab/cufflinks) for cufflinks.

## Cufflinks

[TODO] Brief introduction of cufflinks.

## Installation

### Install Boost

To build Cufflinks, you need to have Boost C++ libraries (version 1.47+) installed. The libraries can be installed at [here](http://www.boost.org/users/download/). 

Download and extract Boost library.

```shell
wget https://dl.bintray.com/boostorg/release/1.65.1/source/boost_1_65_1.tar.gz \
-O boost_1_65_1.tar.gz; tar zxf boost_1_65_1.tar.gz
```

Go to the root directory of Boost library.

```shell
cd boost_1_65_1
```

Run `boostrap.sh` script to make `bjam` executable.

```shell
./bootstrap.sh
```

Run bjam to build Boost. If you are on a 32-bit Linux system, change value of `address_model` argument into 32. 

```shell
bjam --prefix=<YOUR_BOOST_INSTALL_DIRECTORY> \
--toolset=gcc architecture=x86 \
address_model=64 link=static \
runtime-link=static stage install
```

Please remember `<YOUR_BOOST_INSTALL_DIRECTORY>` since it will be used as an argument when building Cufflinks.

### Install Eigen

Also, you need to have [Eigen C++ library](http://eigen.tuxfamily.org/index.php?title=Main_Page) installed. Eigen is a template library for linear algebra in C++.

```shell
wget http://bitbucket.org/eigen/eigen/get/3.3.4.tar.bz2 -O eigen-3.3.4.tar.bz2;
tar xvf eigen-3.3.4.tar.bz2
```

In the directory extracted, there should be `Eigen` directory. Please remember the path to the `Eigen` directory as `<YOUR_EIGEN_INSTALL_DIRECTORY>` since it will be used as an argument when building Cufflinks.

### Install Cufflinks

Download Cufflinks source tarball from [here](http://cole-trapnell-lab.github.io/cufflinks/install/). (You may download the most recent version 2.2.1.)

Extract Cufflinks source tarball and change directory.

```shell
tar zxvf cufflinks-2.2.1.tar.gz; cd cufflinks-2.2.1
```

Befure building, you need to configure your settings.

```shell
./configure --prefix=<YOUR_CUFFLINKS_INSALL_DIRECTORY> \
--with-boost=<YOUR_BOOST_INSTALL_DIRECTORY>
--with-eigen=<YOUR_EIGEN_INSTALL_DIRECTORY>
```

If you successfully configured your settings, you are ready to install Cufflinks.

```shell
make
make install
```

