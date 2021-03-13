Instruction for installing FFTW library:

Download the latest version from [here](http://www.fftw.org/download.html)

```console
$ cd path_to_file
$ sudo ./configure --enable-threads --enable-openmp --with-g77-wrappers
$ sudo make
$ sudo make install
```
