#!/bin/sh
# average multiple afni files together after conversion with mfdftastic
# Monday, December 5, 2005
# mark bolding
# Wednesday, February 8, 2006
# made more generic

echo "usage: $0 basename 02 03 04 06"

prefix=""
postfix=""
fileroot=$1
shift
num=$#

rm calc+orig.HEAD calc+orig.BRIK
image=`ls ${fileroot}${1}+orig.HEAD`
echo $image
3dcalc -a $image -expr " a "
mv calc+orig.HEAD average+orig.HEAD
mv calc+orig.BRIK average+orig.BRIK

shift
for f in $*	
	do
	image=`ls ${fileroot}${f}+orig.HEAD`
	echo $image
	3dcalc -a $image -b average+orig.HEAD -expr " b + a "
	mv calc+orig.HEAD average+orig.HEAD
	mv calc+orig.BRIK average+orig.BRIK
done

3dcalc -nscale -datum short -a average+orig -expr " a / $num "
mv calc+orig.HEAD ${fileroot}average+orig.HEAD
mv calc+orig.BRIK ${fileroot}average+orig.BRIK

echo "$num images averaged." 
