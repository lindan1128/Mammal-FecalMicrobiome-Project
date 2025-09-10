if(@ARGV!=2){
        print "Usage:perl $0 <input.fq> <name list>\n";
        exit;
}

open IN,$ARGV[1] or die "";
while(<IN>){
        chomp;
    
		$list{$_}=1;       
       
        
}
close IN;

open IN,$ARGV[0] or die "";
while(<IN>){
        chomp;
        if(/^\@.+\s/){
                @temp=split;
                $name=$temp[0];
                $name=~s/^\@//;
#                 print "$name\n";
                if(exists$list{$name}){
                	$permi="n";
                }else{
                	$permi="y";
                	print "$_\n";
                }
                $_=<IN>;
                if($permi eq "y"){
                	print ;
                }
                $_=<IN>;
                if($permi eq "y"){
                	print ;
                }
                $_=<IN>;
                if($permi eq "y"){
                	print ;
                }
                
        }
}
close IN;

