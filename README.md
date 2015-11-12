RADAR2NC
========

This is a utility to convert data downloaded from NCDC archive into NetCDF files that iRadar can readily display.

### Requirements ###

On Linux:

* [GCC] GNU C Compiler
* [NetCDF] NetCDF library framework
* [RSL] TRMM Radar Software Library

On Mac:

* [Xcode 6]


[GCC]: http://gcc.gnu.org
[Xcode 6]: https://developer.apple.com/xcode
[NetCDF]: http://www.unidata.ucar.edu/software/netcdf
[RSL]: http://trmm-fc.gsfc.nasa.gov/trmm_gv/software/rsl

NCDC Archive
------------

NOAA (National Oceanic and Atmospheric Administration) provides historical radar data through the NCDC archive accessible through a web portal

http://www.ncdc.noaa.gov/nexradinv/map.jsp


TRMM Radar Software Library
---------------------------

NASA (National Aeronautics and Space Administration) GSFC (Goddard Space Flight Center) provides a software library for ingesting various format of radar data. The library depends on several other libraries so be sure the NetCDF, and HDF5, if you want supported, are installed prior to the installation of the library. The RSL can be downloaded from

http://trmm-fc.gsfc.nasa.gov/trmm_gv/software/rsl/

