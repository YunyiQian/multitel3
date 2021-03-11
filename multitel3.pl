#!/usr/bin/perl
use strict;
# perl code for computing Green's functions using multitel3.
#          Yunyi Qian
#
#
# defaults
my $code = "./multitel3";
my $code2="./rayinformation";
my $dt=0.05;		# sampling rate
my $nt=8192;
my $deg2km = 111.2;		# distances are in km 
my $dir = ".";
my $r0=6371;

my ($model,$stdmodel,$s0_depth,$s_depth,$dist);

# command-line inputs
@ARGV > 1 or die "Usage: multitel3.pl -Mmodel/stdmodel/depth [-D] [-Nnt[/dt]] [-Odir] distances
	-M:  set crust model name, standard earth model and source depth.
	-D:  distances are in degrees (default is km).
	-N:  set time duration ($nt ), dt ($dt s), $nt must be 2^n
	-O:  set output directory name ($dir).\n";

foreach (grep(/^-/,@ARGV)) {
   my $opt = substr($_,1,1);
   my @value = split(/\//,substr($_,2));
   if ($opt eq "D") {
     $deg2km = 1;
   } elsif ($opt eq "N") {
     $nt = $value[0];
     $dt = $value[1] if $#value > 0;
   } elsif ($opt eq "M") {
     $model = $value[0];
     $stdmodel= $value[1];
     $s0_depth = $value[2] if $#value > 0;    #source depth
     $s_depth = $r0*log($r0/($r0-$value[2])) if $#value > 0;  #flat source depth for Haskell
   } elsif ($opt eq "O") {
     $dir = join('/',@value);
   } else {
     printf STDERR "wrong option\n";
     exit(0);
   }
}
my (@ddd) = grep(!/^-/,@ARGV);
my $name1 = "$dir/${model}direct_${s0_depth}";
mkdir($name1,0777) unless -d $name1;
my $name2 = "$dir/${model}core_${s0_depth}";
mkdir($name2,0777) unless -d $name2;

open(RAY,"| $code2");
printf RAY "$s0_depth\n";
printf RAY "0.00\n";       #receivers depth
printf RAY "3.36 5.80\n"   ; #receiver side Vs and Vp value
printf RAY "3.36 5.80 2.6\n";#PP/SS rebound point's Vs, Vp and Rho
printf RAY "$stdmodel\n";
close(RAY);



open(TEL3,"| $code");
print TEL3 "$model\n";
printf TEL3 "%d %f\n",$nt,$dt;
printf TEL3 "%f\n",$s0_depth;
foreach $dist (@ddd) {
  printf TEL3 "%f %f %s\n",0.0,$dist/$deg2km,$dist;
  printf TEL3 "$name1/$dist.grn.\n";
#}
#foreach $dist (@ddd) {
#  printf TEL3 "%f %f %s\n",0.0,$dist/$deg2km,$dist;
  printf TEL3 "$name2/$dist.grn.\n";
}
close(TEL3);
unlink glob "ray*.info";
exit(0);


