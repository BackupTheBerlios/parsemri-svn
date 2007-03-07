#!/usr/bin/perl
print $prefix=$ARGV[0];
print $outputprefix=$ARGV[1];
system "3dcalc -gscale -datum short  -a ${prefix} -expr \'a\' -prefix /tmp/shorttemp";
$info = `3dinfo -short /tmp/shorttemp+orig`;
print "$info";
if ($info =~ /\[\*\s*(\d\.\d+e-\d+)\]/) {
        `rm /tmp/shorttemp*`;
        system "3dcalc -nscale -datum short -a ${prefix} -expr \'a / $1 \' -prefix ${outputprefix} ";
        print $info;
        print "\n";
}
