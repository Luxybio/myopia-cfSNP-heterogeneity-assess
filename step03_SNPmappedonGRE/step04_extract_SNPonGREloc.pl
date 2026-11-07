#!/usr/bin/perl
# Replace <LOCAL_DATA_ROOT>/<SERVER_DATA_ROOT> with your own paths.
%hash=();
$left_dir="<SERVER_DATA_ROOT>/data/Alignment/regulatoryelement-align/EUR-MGIC/EUR-regulatoryelement-MGIC-left.sam";
$right_dir="<SERVER_DATA_ROOT>/data/Alignment/regulatoryelement-align/EUR-MGIC/EUR-regulatoryelement-MGIC-right.sam";
open(I,$left_dir)||die "$!";
while(<I>)
{
chomp($_);
@data=split(/\t/,$_);
$data[1]=~s/256/0/;
$data[1]=~s/272/16/;
$k=$data[0]."\t".$data[1]."\t".$data[2];
$hash{$k}=$data[3];
}

open(IN,$right_dir)||die "$!";
open(O,">><SERVER_DATA_ROOT>/data/Alignment/regulatoryelement-align/EUR-MGIC/EUR-regulatoryelement-MGIC-relocate.txt")||die "$!";

while(<IN>)
{
 chomp($_);
 @lnc=split(/\t/,$_);
 $lnc[1]=~s/256/0/;
 $lnc[1]=~s/272/16/;
 $sam=$lnc[0]."\t".$lnc[1]."\t".$lnc[2];
 while(my ($key,$value)=each %hash)
 {
  if($sam eq $key)
  {
   $cha=$lnc[3]-$hash{$key};
   if(abs($cha)==31)
     {
    print O $sam."\t".$lnc[3]."\t".$hash{$key}."\t".$cha."\n";
     }
  }
 }
}
close(I);
close(IN);
close(O);
