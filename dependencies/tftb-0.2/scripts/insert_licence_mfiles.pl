#!/usr/bin/perl

#-------------------------------------------------
#the script changes the confidential text by 
#that of the license CeCILL .
#
#-------------------------------------------------

use strict;

my $tftb_root = "mfiles";
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
        if ($_ =~ /\.m$/i) 
	{ 	
		my $mo="CONFIDENTIAL PROGRAM";
		my $nblinet = 1;
		my $nbline = 1;
		my $line;
		
		open (FILE, "<$tftb_root/$_") || die "impossible to open $_.";	
		while($line=<FILE>)
		{
			$nblinet++;
		}
		close (FILE);
		
		
		open (TEMP, ">>$tftb_root/temp.log") || die "impossible to open temp.log";
		open (FILE, "<$tftb_root/$_") || die "impossible to open $_.";
		while(($line=<FILE>) and ($line !~ /$mo/))
		{
			print TEMP $line;
			$nbline++;
		}
		close (FILE);
	
	
		open (ENTETE, "< licence") || die "impossible to open entete";
		while (my $line = <ENTETE>)
	   	{
			print TEMP "$line";
	   	}
	     	close (ENTETE);
		
		
		$nbline=$nbline+3;
		open (FILE, "<$tftb_root/$_") || die "impossible to open $_.";
		my @tab=<FILE>;
		for($nbline..$nblinet)
            	{
			print TEMP "$tab[$_]";	
		}	
	  	
		close (FILE);
		close (TEMP);
		
		
		`rm $tftb_root/$_`;
		`mv $tftb_root/temp.log $tftb_root/$_`;					
	}
    }
}
