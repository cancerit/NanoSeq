gzstream
========

This is only a mirror of the original gzstream-library developed by The Computational Geometry Group at UNC Chapel Hill. For more information about this lib, see the [cs.unc.ed homepage](http://www.cs.unc.edu/Research/compgeom/gzstream/)

installation
------------

to install please make sure you have zlib installed. In ubuntu its included in package `zlib1g-dev`.  
Download the project or clone it with
```
git clone https://github.com/kanedo/gzstream.git
```
  
change into the dir
```
cd gzstream
```
now type 
```
sudo make install
``` 
to install gzstream.

This will will compile the gzstream library, test it and install it to the folllowing locations:  
`/usr/include/gzstream.h`  
`/usr/lib/libgzstream.a`  

remove library
--------------

To remove the library just use the make file and type
```
sudo make uninstall
```

usage
-----

To use the library with `Make` make sure you add `-lgzstream` and `-lz` to your `LDFLAGS`.

FAQ
---
I get `undefined reference to gzstream` error:  
make sure you add the LDFlags at the end
