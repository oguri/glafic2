#!/usr/bin/perl
#
# convert galfile for glafic to ds9 region file
#

# size of circle in arcsec (for WCS)
$rr=0.8;
# color of cirle
$clr='blue';

if((($#ARGV+1)!=1)&&(($#ARGV+1)!=2)){
    die "usage:\n > conv_galfile_dat2reg.pl model_hoge.input [galfile_hoge.dat]\n    [galfile_hoge.dat] is optional\n";
}

## check model_hoge.input file
open(IN,"$ARGV[0]") || die "File Not Found\n";

$x0="";
$y0="";
while(<IN>){
    if(substr($_,0,1)  ne '#'){
        @data=split;
        
        $com=$data[0];
	if($com eq 'prefix'){ $prefix=$data[1]; }
        if($com eq 'wcs_ra0'){  $x0=$data[1]; }
        if($com eq 'wcs_dec0'){ $y0=$data[1]; }
	if($com eq 'galfile'){ $ofile=$data[1]; }
    }
}

close(IN);

if($x0 =~ /[0-9]/){
    printf "use WCS (ra0: %s , dec0 %s)\n",$x0,$y0; 
} else {
    printf "use relative coordinates\n";
}

## read galfile, write reg file
if(($#ARGV+1)==1){
    $fc=$ofile;
} else {
    $fc=$ARGV[1];
}

open(IN,"$fc") || die "File Not Found\n";

$f1="reg_".substr($fc, rindex($fc, '/') + 1);;

$f1 =~ s/.dat/.reg/g;

printf("input galfile     : %s\n",$fc);
printf("output region file: %s\n",$f1);

open(OUT,">$f1") || die "Error: File Not Found\n";

print OUT "# Region file format: DS9 version 4.1\n";
print OUT "global color=green font=\"helvetica 8 normal roman\"\n";
if($x0 =~ /[0-9]/){
    print OUT "fk5\n";
} else {
    print OUT "physical\n";
}

open(IN,"$fc") || die "File Not Found\n";

while(<IN>){
    if(substr($_,0,1)  ne '#'){
	@data=split;
	$dr=$data[0];
	$dd=$data[1];
	
	if($x0 =~ /[0-9]/){
	    $dec=(3.141592/180.0)*$y0;
	    $ys=$y0+$dd/3600.0;
	    $xs=$x0-$dr/(3600.0*cos($dec));
	} else {
	    $xs=$dr;
	    $ys=$dd;
	}
		    
	if($x0 =~ /[0-9]/){
	    print OUT "circle($xs,$ys,$rr\") # width=2 text={} color=$clr\n";
	} else {
	    print OUT "circle($xs,$ys,$rr) # width=2 text={} color=$clr\n";
	}
    }
}

close(IN);
close(OUT);

