function [ cvx_ver2, cvx_bld ] = cvx_version
cvx_ver = 1.21;
cvx_bld = '795';
cvx_bdate = '2010-05-31 08:26:19';
cvx_dbld = '795';
cvx_ddate = '2010-05-31 08:26:19';
if nargout == 0,
   fprintf( '\n' );
   fprintf( 'CVX version %g\n', cvx_ver );
   fprintf( '    Code: build %s, %s\n', cvx_bld, cvx_bdate );
   fprintf( '    Documentation: build %s, %s\n', cvx_dbld, cvx_ddate );
   if exist( 'OCTAVE_VERSION', 'var' ),
       fprintf( 'GNU Octave %s on %s\n', version, computer );
       fprintf( 'NOTE: Sorry, Octave support is not yet functional.\n' );
   else
       verd = ver('MATLAB');
       fprintf( 'MATLAB version %s %s on %s\n', verd.Version, verd.Release, computer );
   end
   fprintf( '\n' );
else
    cvx_ver2 = cvx_ver;
end

