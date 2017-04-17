#!/usr/bin/perl
use warnings;

# Print a list  of all subdirectories into an array via locate
my @Directories = `locate *.xtc | grep APO-data`;

# make a txtfile to put the modified directory paths in
my $newfile="list.txt";
open (LIST, '>', "$newfile") || die $!;

# remove unwanted parts to allow for uniq "sorting"
foreach $Directory (@Directories){
        my @Folders = split(/\/0/, "$Directory"); # split the filepath at each /0 (not /00 because 010)
        print LIST "$Folders[0]"; # print this new array to the txtfile
        print LIST "\n"; # put newlines in between each array printout

}
close (LIST); # save the file?

# make new list by removing repeat directories
`less list.txt | uniq > uniq.txt`;
open (NEWLIST, '<', "uniq.txt");
chomp (my @UniqDirectories = <NEWLIST>);
close (NEWLIST);

print scalar (@Directories), " total .xtc paths in List 1.\n"; # should equal how many total .tpr paths there are
print scalar @UniqDirectories, " unique .xtc paths in List 2.\n";  # how many uniq .tpr types are there?


exit;
