This directory contains MKL extensions for ALGLIB.

ALGLIB may utilize Intel MKL to speed-up some  calculations.  It  may  use
already existing MKL installation  or special lightweight MKL distribution
which is shipped with  Commercial Edition  of  ALGLIB  and  includes  only
those MKL functions which are used by ALGLIB.

IMPORTANT: you do not have to download/install Intel MKL in order  to  use
           it with ALGLIB. All the necessary code is already included into
           files in this directory. Connecting them to ALGLIB is  easy and
           fast!

This directory contains just precompiled MKL extensions; you still have to
compile the rest of ALGLIB yourself. Sections below contain information on
compiling/linking under Windows and Linux.  You  can  find  more  detailed
discussion  in  the  'Linking with Intel MKL' section of ALGLIB  Reference
Manual.

NOTE: if you activate MKL extensions, then  following  license  terms  are
      applied to you (in addition to ALGLIB License Agreement):
      * mkl-2021-license.txt
      * mkl-2021-third-party-programs.txt



=== USING MKL EXTENSIONS ON 64-BIT WINDOWS ===============================

1. compile ALGLIB source files with following preprocessor symbols defined
   at the global level:
   * AE_MKL                 to use Intel MKL functions in ALGLIB
   * AE_OS=AE_WINDOWS       [optional] to activate multicore capabilities
   * AE_CPU=AE_INTEL        [optional] to tell AGLIB that it is Intel/AMD

2. link your application with ALGLIB and with following import library:
   * alglibXYZ_64mkl.lib    (this LIB is an import library for larger DLL)

3. place following DLL file  into directory where it can  be  found during
   application startup (usually - application dir):
   * alglibXYZ_64mkl.dll    (this DLL contains MKL code)



=== USING MKL EXTENSIONS ON 64-BIT LINUX =================================

1. compile ALGLIB source files with following preprocessor symbols defined
   at the global level:
   * AE_MKL                 to use Intel MKL functions in ALGLIB
   * AE_OS=AE_POSIX         [optional] to activate multicore capabilities
   * AE_CPU=AE_INTEL        [optional] to tell AGLIB that it is Intel/AMD

2. link your application with ALGLIB and with following shared library:
   * libalglibXYZ_64mkl.so  (this SO contains MKL code)
   
3. make sure that shared library can be found at runtime, i.e. either:
   (a) place libalglibXYZ_64mkl.so into directory where it  can  be  found
       during binary startup (/lib, /usr/local/lib, and so on)
   (b) change LD_LIBRARY_PATH environment variable to include app dir
   (c) specify rpath during application compilation



=== USING MKL EXTENSIONS ON 32-BIT WINDOWS ===============================

Exactly  same  as  using  on  64-bit  Windows,  except  for  the fact that
alglibXYZ_64mkl.lib/dll becomes alglibXYZ_32mkl.lib/dll



=== USING ALGLIB WITH YOUR OWN INSTALLATION OF INTEL MKL =================

If you already have Intel MKL installed on your system and want ALGLIB  to
be  compiled  against  your  own  installation,  then  you  should perform
following steps:

1. make  sure  that  your  installation  of  Intel MKL can be seen by your
   compiler (all environment variables are properly configured)  and  your
   link line contains all the necessary parameters. We  recommend  you  to
   use Intel MKL Link Line Advisor tool from Intel website.  We  recommend
   you to use ILP64 MKL interface on 64-bit systems, because it has better
   interoperability with ALGLIB, although everything will  work  with  any
   kind of interface.

2. compile mkl4alglib.c with the rest of ALGLIB, as described in step (3).
   This file serves as interface layer between ALGLIB and Intel MKL.

3. compile ALGLIB source files  (including  mkl4alglib.c)  with  following
   preprocessor symbols defined at the global level:
   * AE_MKL                 to use Intel MKL functions in ALGLIB
   * AE_OS=AE_POSIX         [mandatory, on Linux]
   * AE_OS=AE_WINDOWS       [mandatory, on Windows]
   * AE_CPU=AE_INTEL        [optional] to tell AGLIB that it is Intel/AMD
   mandatory OS definitions are necessary  for  ALGLIB-MKL  interface  to
   know about OS it is running on.
