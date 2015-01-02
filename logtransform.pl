use strict;

#logtransform.pl combined_nonred.txt combined_log.txt

my $usage = "$0 input output [-nolog]\n";
my $arg=0;
my $input_file = $ARGV[$arg++] or die $usage;
my $output_file = $ARGV[$arg++] or die $usage;

open (OUTPUT, ">$output_file") or die $!;
open (INPUT, $input_file) or die $!;
my $line = <INPUT>;
my @currentline = split(/\t/, $line);
for my $i (0..$#currentline)
{
    if ($currentline[$i] =~ m/Reference/)
    {
        next;
    }
    else
    {
        print OUTPUT "$currentline[$i]\t";
    }
}
print OUTPUT "\n";
$line =~ s/\n$//;
my $log10 = log(10);
my $count=0;
while(my $line = <INPUT>){
    $line =~ s/\n$//;
    my($oligo1,@data1)=split(/\t/, $line);
    #if($oligo1 ne "A_23_P1473"){ next; }
    my @grn1=();
    my @red1=();
    my ($mr,$mg)=(0,0);
    for(my $i=0; $i<@data1; $i+=2){
        my $z = log($data1[$i]);
        push(@grn1,$z);
        $mg += $z;
        $z = log($data1[$i+1]);
        push(@red1,$z);
        $mr += $z;
    }
    $mr /= @red1;
    $mg /= @red1;
    print OUTPUT $oligo1;
    my ($r,$vx,$vy) = correlation(\@grn1,\@red1);
    my $slope = $r*sqrt($vy/$vx);
    if($mr-$mg < -1.1 && $slope > 0.3){    #1.1=ln(3)
        for(my $i=0; $i<@grn1; $i++){
            my $x = int(10000*$grn1[$i]/$log10+0.5)/10000;
            print OUTPUT "\t$x";
        }
    }else{
        for(my $i=0; $i<@grn1; $i++){
            my $x = int(10000*($grn1[$i]-$red1[$i]+$mr)/$log10+0.5)/10000;
            print OUTPUT "\t$x";
        }
    }
    print OUTPUT "\n";
}
close INPUT;
close OUTPUT;
exit(0);


#**************************
sub   correlation
#**************************
{
my $xr = shift;
my $yr = shift;
my $n = shift;

my $sx = 0;
my $sxx = 0;
my $sxy = 0;
my $sy = 0;
my $syy = 0;
if(!$n){
    $n = @$xr;
}
for(my $i=0; $i<$n; ++$i){
    my $x = $xr->[$i];
    my $y = $yr->[$i];
    $sx += $x;
    $sxx += $x*$x;
    $sxy += $x*$y;
    $sy += $y;
    $syy += $y*$y;
}
my $vx = $sxx-$sx*$sx/$n;
my $vy = $syy-$sy*$sy/$n;
my $r = $sxy-$sx*$sy/$n;
if($vx*$vy >0){
    $r = $r/sqrt($vx*$vy);
}
$vx /= ($n-1);
$vy /= ($n-1);
return ($r,$vx,$vy);
}

#**************************
sub   median
#**************************
{
my $xr = shift;
my $n = @$xr;
if(!$n){ return 0; }
my $median = 0;
my @sorted = sort {$a<=>$b} @$xr;
my $i = $n/2;
if($i > int($i)){
    $i = int($i);
    $median = $sorted[$i];
}else{
    $median = ($sorted[$i] + $sorted[$i-1])/2;
}
return $median;
}

#**************************
sub   covariance
#**************************
{
my $xr = shift;
my $yr = shift;

my $sx = 0;
my $sxy = 0;
my $sy = 0;
my $n = @$xr;
for(my $i=0; $i<$n; ++$i){
    my $x = $xr->[$i];
    my $y = $yr->[$i];
    $sx += $x;
    $sxy += $x*$y;
    $sy += $y;
}
return $sxy-$sx*$sy/$n;
}

#***********************************
sub regression_oneway
#***********************************
{
my $y_ref = shift;
my $x_ref = shift;

my $sx=0;
my $sy=0;
my $syx=0;
my $sxx=0;
my $syy=0;
my $n=0;
for(my $i=0; $i<@$x_ref; ++$i){
    $sx += $x_ref->[$i];
    $sy += $y_ref->[$i];
    $n++;
}
if($n < 2){
    print "No pattern defined\n";
    exit(0);
}
$sx /= $n;
$sy /= $n;
for(my $i=0; $i<@$x_ref; ++$i){
    my $x = $x_ref->[$i]-$sx;
    my $y = $y_ref->[$i]-$sy;
    $syx += $x*$y;
    $sxx += $x*$x;
    $syy += $y*$y;
}
my $b1 = $syx/$sxx;
my $a1 = $sy - $b1*$sx;
return ($a1,$b1,$sy,$sx,$syy/($n-1),$sxx/($n-1));
}


