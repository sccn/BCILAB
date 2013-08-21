function osc_make()

  if isunix

    % requires liblo to be available on your system, e.g. in /usr/local/lib or /usr/lib.

    cd(fileparts(which('osc_make')))

    mex -llo osc_new_address.c
    mex -llo osc_free_address.c
    mex -llo osc_new_server.c
    mex -llo osc_free_server.c
    mex -llo osc_send.c
    mex -llo osc_recv.c

  else

    % compile commands to build the osc library DLLs

    % these dlls are built with gcc / MinGW
    % gnumex is used to enable gcc support for mex.
    % the liblo object files are also built with gcc.
    % obviously you will need to change the paths etc for your system.

    cd y:\matlab\win\osc\

    mex -Iy:\matlab\win\osc\liblo-0.22 -f y:\matlab\win\gnumex\mexopts.bat y:\matlab\osc\osc_new_address.c y:\matlab\win\osc\liblo-0.22\src\.libs\*.o c:\mingw\lib\libws2_32.a
    mex -Iy:\matlab\win\osc\liblo-0.22 -f y:\matlab\win\gnumex\mexopts.bat y:\matlab\osc\osc_free_address.c y:\matlab\win\osc\liblo-0.22\src\.libs\*.o c:\mingw\lib\libws2_32.a
    mex -Iy:\matlab\win\osc\liblo-0.22 -f y:\matlab\win\gnumex\mexopts.bat y:\matlab\osc\osc_new_server.c y:\matlab\win\osc\liblo-0.22\src\.libs\*.o c:\mingw\lib\libws2_32.a
    mex -Iy:\matlab\win\osc\liblo-0.22 -f y:\matlab\win\gnumex\mexopts.bat y:\matlab\osc\osc_free_server.c y:\matlab\win\osc\liblo-0.22\src\.libs\*.o c:\mingw\lib\libws2_32.a
    mex -Iy:\matlab\win\osc\liblo-0.22 -f y:\matlab\win\gnumex\mexopts.bat y:\matlab\osc\osc_send.c y:\matlab\win\osc\liblo-0.22\src\.libs\*.o c:\mingw\lib\libws2_32.a
    mex -Iy:\matlab\win\osc\liblo-0.22 -f y:\matlab\win\gnumex\mexopts.bat y:\matlab\osc\osc_recv.c y:\matlab\win\osc\liblo-0.22\src\.libs\*.o c:\mingw\lib\libws2_32.a
    
    cd ..\..\osc
    
  end
