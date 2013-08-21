#!/usr/bin/perl

#-------------------------------------------------
#
# The script execute the tests under octave and 
# sent the errors in error.log
#
#-------------------------------------------------

use strict;

my $tftb_root = "../tests";
@ARGV = "$tftb_root";

my $dir = $ARGV[0];

ListRep($dir);

sub ListRep {
     my ($dir) = @_;
     if (! -e $dir ) {
     print "Unknown directory ($dir).";
     return undef;
     }
   
     if (! -d $dir ) {
     print "$dir is not a directory";
     return undef;
     }
        
     if (! opendir( DIR, $dir) ) {
     print "Impossible to open the directory $dir : $!.";
     return undef;
     }
     
     my @files = grep !/(?:^\.$)|(?:^\.\.$)/, readdir DIR;
     closedir DIR;
     
     my $date =localtime(time);
     my $version = `octave -version`;     
     open (FILE, ">>error.log") || die "impossible to open error.log";
     print FILE "Date the file:$date\n";
     print FILE "\n";
     print FILE "$version\n";
     close (FILE);
      
     print "\nFiles contained\n";
     foreach(@files) {
        print $_."\n";
        if ($_ =~ /\.m$/i) 
	{ 
		unless (fork)
		{
		
		open (FILE, ">>error.log") || die "impossible to open error.log";
		print FILE "\n";
		print FILE "NAME OF THE FILE : $_\n";
        	close (FILE); 
		
		system( "octave -q -p ../mfiles $tftb_root/$_  2>> error.log");
		 	
		}
		else {
  		die "\n";
		}
	}
    }
}
 
