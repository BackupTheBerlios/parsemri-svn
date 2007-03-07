#!/usr/bin/perl -w 
# November, 2005
# mark bolding

$,=" ";
$\="\n";
$debug = 1;

if ( scalar(@ARGV) < 1 ) { 
	print "usage: fdftastic2.pl shots prefix folder \n shots is the number of shots per volume \n prefix is prefixed to output files \n folder is the .img folder with the .fdf files";
	#$numargs = scalar(@ARGV);
	#print $numargs; 
	#print $ARGV[0];
	die "arrrr... \n";
}
#---------------use line args.
$shots = $ARGV[0]; #number of shots not in fdf, maybe in procpar?

$outputprefix = "";
if ( $ARGV[1] ) {
	$outputprefix = $ARGV[1];
}
else {
	print "using defalt filename...";
	$outputprefix = "functional_volume.nii.gz";
}

if ( $ARGV[2] ) {
	chdir $ARGV[2];  # i'm goin in, cover me...
}
else {
	print "working in current directory... \n";
}

#--------------- get files
@slicefiles = glob "*.fdf";
$numslicefiles = scalar(@slicefiles);

if ($slicefiles[0] =~ /(\w+)\.(fdf)$/) {
    $rootname=$1;
	$extension=$2;
}

open( FILE, "< $slicefiles[0]" ) or die "Can't open $slicefiles[0] : $!";
		
# ------------parse first .img header
while (<FILE>) {
    chomp;	# strip record separator
    
    # ROI
    if (/float\s*roi\[\] = \{(\d+\.\d+),(\d+\.\d+),(\d+\.\d+)\}/) {
    	$ro = $1 * 10;
    	$pe = $2 * 10;
    	$thk = $3 * 10;
    }
    
    if (/float  TR = (\d+\.\d+)/) {
    	$tr = $1 * $shots; #tr in fdf file is per shot, erk.
	}
	
	if (/int    slices = (\d+)/) {
		$slicespervolume = $1;
	
	}
	
	if (/float  matrix\[\] = \{(\d+), (\d+)\}/) {
		$matrixX = $1;
		$matrixY = $2;
	}
	
	if (/float  orientation\[\] = \{(.*)\}/) {
		if ( $1 eq "-1.0000,0.0000,-0.0000,0.0000,-1.0000,0.0000,-0.0000,0.0000,1.0000" ) {$orient = "axial"}
		else {$orient = "Unknown"}
	
	}
	
	if (/float  location\[\] = \{(-?\d+\.\d+),(-?\d+\.\d+),(-?\d+\.\d+)\}/) {
		$firstsliceposition = $3 * 10;
	
	}
	
	
}

close FILE;

#open last file
open( FILE, "< $slicefiles[$numslicefiles - 1]" ) or die "Can't open $slicefiles[$numslicefiles - 1] : $!";
		
# ------------parse first .img header
while (<FILE>) {
    chomp;	# strip record separator

	
	if (/float  location\[\] = \{(-?\d+\.\d+),(-?\d+\.\d+),(-?\d+\.\d+)\}/) {
		$lastsliceposition = $3 * 10;
	
	}
	
	
}

close FILE;

# -----------calculate values and args for AFNI to3d to convert fdf files
#change to use first and last slice position moved out by thk/2
#A is -P, L is -R etc. up is -down, good is -bad
#new vars x1 x2 y1 y2 z1 z2  instead of 
$zslab = $thk * $slicespervolume;
$sliceoffset = $firstsliceposition - ($thk / 2) + ($zslab / 2);

#calculate FOV args for axial orientation, ignore readout offset
$ro = $ro/2; $pe = $pe/2;  $inferior = $zslab/2 + $sliceoffset; $superior = $zslab/2 - $sliceoffset;

$x1 = "${pe}L";
$x2 = "${pe}R";
$y1 = "${ro}A";
$y2	= "${ro}P";
$z1 = "${superior}S";
$z2 = "${inferior}I";

$FOVargs = "-xFOV ${x1}-${x2} -yFOV ${y1}-${y2} -zFOV ${z1}-${z2}";

# name of output file
$prefix = "output.nii.gz";
$prefixargs = "-prefix $prefix";

# total number of 2d images 
$totalslices = $numslicefiles;
$volumes = $totalslices / $slicespervolume;

# generate time args for func or anat
$timeargs = "";
if ($volumes > 1) {
	print "$volumes timepoints";
	$timeargs = "-time:tz $volumes $slicespervolume ${tr}ms zero";
#	$timeargs = "-time:zt $slicespervolume $volumes ${tr}ms zero";
}

#file specification
$filespec = "3Df:-1:0:$matrixX:$matrixY:1:*.$extension";

if ($debug == 1) {
	print "-------------\nmatrix = $matrixX by $matrixY \n";
	print "first file = $rootname.$extension \n";
	#print "@slicefiles \n";
	print " number of files = $numslicefiles";
	print "  tr = $tr";
	print "  readout = ${ro} mm , phase encode = ${pe} mm , thk = ${thk} mm \n";
	print "slices per volume = $slicespervolume \n";
	print "orientation = $orient \n ";
	print "first slice position = $firstsliceposition mm , last slice position = $lastsliceposition , zslab thickness = $zslab mm , slice offset = $sliceoffset mm \n";

}


# to3d command----------------
$command="to3d $prefixargs $timeargs $FOVargs $filespec";
print "$command";
`$command > to3dout.txt 2>&1`; # do it!

#find scaling factor for short conversion. why doesnt -datum short work?
system "3dcalc -gscale -datum short  -a ${prefix} -expr \'a\' -prefix shorttemp.nii.gz"; 
$info = `3dinfo -short shorttemp.nii.gz`;
print "$info";
if ($info =~ /\[\*\s*(\d\.\d+e-\d+)\]/) {
	`rm shorttemp.nii.gz`;
	system "3dcalc -nscale -datum short -a ${prefix} -expr \'a / $1 \' -prefix ${outputprefix} ";
	`rm ${prefix}`;
	print $info;
	print "\n";
}
	