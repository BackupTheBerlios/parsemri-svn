#!/usr/bin/perl -w
#3dPARtoANALYZE  somefile.par 


$,=" ";
$\="\n";
$debug = 1;
 
if ($ARGV[0] =~ /([A-Za-z_0-9]+)\.([Pp][Aa][Rr])$/) {
        $rootname=$1;
	$extension=$2;
}
		
# parse PAR file.
while (<>) {
    chomp;	# strip record separator
    # field of view
    if (/FOV \((..),(..),(..)\) \[.*\].*:\ *([0-9]+\.[0-9]+)\ *([0-9]+\.[0-9]+)\ *([0-9]+\.[0-9]+)/) {
	$fov{$1}=$4;
	$fov{$2}=$5;
	$fov{$3}=$6;
    }
    # matrix size
    if (/Recon resolution \((.), (.)\).*:\ *([0-9]+)\ *([0-9]+)/) {
	$reconres{$1}=$3;
	$reconres{$2}=$4;
    }
    # angulation
    if (/Angulation midslice\((..),(..),(..)\)\[.*\].*:\ *(-?[0-9]+\.[0-9]+)\ *(-?[0-9]+\.[0-9]+)\ *(-?[0-9]+\.[0-9]+)/) {
	$angulation{$1}=$4;
	$angulation{$2}=$5;
	$angulation{$3}=$6;
    }
    # distance from origin of FOV center
    if (/Off Centre midslice\((..),(..),(..)\) \[.*\].*:\ *(-?[0-9]+\.[0-9]+)\ *(-?[0-9]+\.[0-9]+)\ *(-?[0-9]+\.[0-9]+)/) {
	$origin{$1}=$4;
	$origin{$2}=$5;
	$origin{$3}=$6;
    }
    # slice thickness
    if (/Slice thickness.*:\ *([0-9]+\.[0-9]+)/) {
    	$thk=$1;
    }
    # gap
    if (/Slice gap.*:\ *([0-9]+\.[0-9]+)/) {
    	$gap=$1;
    }
    # number of slices
    if (/number of slices.*:\ *([0-9]+)/) {
	$slices=$1;
    }
    # image bit depth
    if (/Image pixel size.*:\ *([0-9]+)/) {
	$bitdepth=$1;
    }	
    # TR
    if (/Repetition time.*:\ *([0-9]+\.[0-9]+)/) {
    	$tr=$1;
    }
    # number of volumes
    if (/number of dynamics.*:\ *([0-9]+)/) {
	$volumes=$1;
    }	
 
}  

if ($debug == 1) {
	print "$rootname.$extension";
	print "--FOV"; 
	print "ap=$fov{ap}";
	print "fh=$fov{fh}";
	print "rl=$fov{rl}";
	print "--recon matrix";
	print "x=$reconres{x}";
	print "y=$reconres{y}";
	print "--angulation";
	print "ap=$angulation{ap}";
	print "fh=$angulation{fh}";
	print "rl=$angulation{rl}";
	print "--origin";
	print "ap=$origin{ap}";
	print "fh=$origin{fh}";
	print "rl=$origin{rl}";
    	print "--slice thickness = $thk";
    	print "--slice gap = $gap";
	print "--slices = $slices";
    	print "--bit depth = $bitdepth";
	print "--TR = $tr";
}

# -----------calculate values and args for AFNI to3d to convert REC file
# 
$prefixargs = "-prefix $rootname";
# generate time args for func or anat
$timeargs = "";
$totalslices = $slices;
$avwargs = "-4D -orient LPI";
if ($volumes > 1) {
	print "$volumes timepoints";
	$timeargs = "-time:tz $volumes $slices ${tr}ms zero";
#	$timeargs = "-time:zt $slices $volumes ${tr}ms zero";
	$totalslices=$volumes*$slices;
}
#calculate FOV args
$anterior = abs($origin{ap}-$fov{ap}/2);
$posterior= abs($origin{ap}+$fov{ap}/2);
$foot     = abs($origin{fh}-$fov{fh}/2);
$head     = abs($origin{fh}+$fov{fh}/2);
$right    = abs($origin{rl}-$fov{rl}/2);
$left     = abs($origin{rl}+$fov{rl}/2);

$FOVargs = "-xFOV ${right}R-${left}L -yFOV ${anterior}A-${posterior}P -zFOV ${foot}I-${head}S";

#file specification
$recextension="rec";
if ($extension =~ /PAR/) {
	$recextension="REC";
}
$filespec = "3D:-1:0:$reconres{x}:$reconres{y}:$totalslices:$rootname.$recextension";

# to3d command
print $command="to3d $prefixargs $timeargs $FOVargs $filespec";
system $command; 
# convert to analyze command
print $command="3dAFNItoANALYZE $avwargs $rootname ${rootname}+orig";
system $command;
