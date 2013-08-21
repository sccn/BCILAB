#!/usr/bin/perl

use strict;

my $tftb_root = "refguide";
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
      
     print "\nFiles contained\n";
     foreach(@files) {
        print $_."\n";
        if ($_ =~ /\.tex$/i) 
	{ 	
		my $line;
				
		open (TEMP, ">>$tftb_root/temp.log") || die "impossible to open temp.log";
		open (ENTETE, "< licence_doc") || die "impossible to open licence_doc";
		while (my $line = <ENTETE>)
	   	{
			print TEMP "$line";
	   	}
	     	close (ENTETE);
		
		
		open (FILE, "<$tftb_root/$_") || die "impossible to open $_.";
		while (my $line = <FILE>)
            	{
			print TEMP "$line";	
		}	
	  	
		close (FILE);

		close (TEMP);
		
		
		`rm $tftb_root/$_`;
		`mv $tftb_root/temp.log $tftb_root/$_`;					
	}
    }
}
