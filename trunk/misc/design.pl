#!/usr/bin/perl
# Monday, December 12, 2005 5:20:55 PM
# mark bolding
# spit out an EV vector
print "usage: design delay on off total\n" if ($#ARGV < 3);
$delay=$ARGV[0];
$on   =$ARGV[1];
$off  =$ARGV[2];
$total=$ARGV[3];
#print "delay:$delay on:$on off:$off total:$total\n";
$count = 0;
$i = 0;
while ($i < $delay && $count < $total) 
{
	print "0 ";
	$count++; #print "$count ";
	$i++;
}
print"\n"; 
while ($count < $total)
{
	$i = 0;
	while ($i < $on && $count < $total) 
	{
		print "1 ";
		$count++; #print "$count ";
		$i++;
	}
	print"\n"; 
	$i = 0;
	while ($i < $on && $count < $total) 
	{
		print "0 ";
		$count++; #print "$count ";
		$i++;
	}
	print"\n"; 
}