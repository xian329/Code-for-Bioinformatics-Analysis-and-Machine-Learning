

use strict;
#use warnings;

use File::Copy;

my $newDir="files";
unless(-d $newDir){
	mkdir $newDir or die $!;
}

opendir(RD, ".") or die $!;
my @allFiles=readdir(RD);
closedir(RD);

foreach my $subDir(@allFiles)
{	
	next if($subDir eq '.');
	next if($subDir eq '..');
	if((-d $subDir) && ($subDir ne $newDir))
	{
		opendir(SUB,"./$subDir") or die $!;
		while(my $file=readdir(SUB))
		{
			if($file=~/\.tsv$/)
			{
				#`cp ./$subDir/$file ./$newDir`;
				copy("$subDir/$file","$newDir") or die "Copy failed: $!";
			}
		}
		close(SUB);
	}
}

