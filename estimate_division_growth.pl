#!/import/bc2/soft/bin/perl5/perl

$mintopthresh = 100; # position cutoff at channel exit

$infile = shift(@ARGV);
open(F,$infile) || die "cannot open infile $infile\n";
#get first comment line
$_ = <F>;
my @time = ();
my @size = ();
my @florsize = ();
my @fluo = ();
my $t = 0;


#Space for residuals
my @resbegin = ();
my @resend = ();
my @countbegin = ();
my @countend = ();
my @resbeginbot = ();
my @resendbot = ();
my @countbeginbot = ();
my @countendbot = ();
for($i=0;$i<100;++$i){
    $countbegin[$i] = 0;
    $countend[$i] = 0;
    $resbegin[$i] = 0;
    $resend[$i] = 0;

    $countbeginbot[$i] = 0;
    $countendbot[$i] = 0;
    $resbeginbot[$i] = 0;
    $resendbot[$i] = 0;
}


print ">GROWTH\n";
print "#laneID\tcellID\tparID\tdaughter-type\tstart-time\tend-time\tend-type\tismintoppixel\tavpos\tavnumfrombottom\tlast-time-fit\tRsize_vs_time\tslope_size_vs_time\terror-slope\tlogsizestart\tlogsize_end\tlogsize-errobar\tdlogsize\tsigmadlogsize\tmeasfluostart\tmeasfluoend\tRfluo_vs_size\tfluo_per_micron\tfluo_zero_micron\testfluostart_size\testfluoend_size\tRfluo_vs_time\tdfluo_per_sec\testfluostart_time\testfluoendtime\n";

#Only use cells that do not go below pixel 100!!

while(<F>){
    #new cell
    if($_ =~ /\>CELL/){
	$line = $_;
	#print STDERR "doing $line";
	#print STDERR "previous has time $t endtype $endtype and mintop $mintop\n";
	#check if there is a previous example to parse
	if($t<4 || $mintop < $mintopthresh){ # || $endtype ne "div"
		if($t>0){
	    print "$celline NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA\n";
		}
	} else {
# 	if($t>=4 && $mintop >= $mintopthresh){
	    #number of time points
	    $T = @time;
	    
	    #at least take 2/3 of the point
	    $Tmin = int($T * 2/3);
	    #print "Tmin is $Tmin\n";
	    $Lmax = -10000000000.0;
	    for($Tend=$Tmin;$Tend<=$T;++$Tend){
		#get the linear fit
		$musize = 0;
		$mut = 0;
		$vart = 0;
		$varsize = 0;
		$covar = 0;
		for($t=0;$t<$Tend;++$t){
		    $curtime = $time[$t]-$time[0];
		    $mut += $curtime;
		    $musize += $size[$t];
		    $vart += $curtime*$curtime;
		    $varsize += $size[$t]*$size[$t];
		    $covar += $curtime*$size[$t];
		}
		
		if($Tend <= 0){
		    print STDERR "PROBLEM: Tend $Tend at t $t, T is $t, Tmin $Tmin\n"
		}
		
		
		$mut /= $Tend;
		$musize /= $Tend;
		$varsize /= $Tend;
		$covar /= $Tend;
		$vart /= $Tend;
		$covar  -= $mut*$musize;
		$vart -= $mut*$mut;
		$varsize -= $musize*$musize;

		$alpha = $covar/$vart;
		#Check for zero variance in size
		if($varsize <= 0){
		    $varsize = 0.5;
		}
		
		$r = $covar/sqrt($vart*$varsize);
		$size0 = $musize - $alpha * $mut;

		#Get variance
		$sigsq = 0;
		for($t=0;$t<$Tend;++$t){
			$d = $size[$t]-$size0-$alpha*($time[$t]-$time[0]);
			$sigsq += $d*$d;
		}
		$D = $sigsq;
		$sigsq /= $Tend;

		#now get likelihood
		if($sigsq <=0){
			$sigsq = 0.01;
		}

		$L= -0.5*$Tend*log($sigsq);
		$L -= 0.5*$D/$sigsq;#from the variance
		#from the endpoint
		$L -= ($T-$Tend)*log($maxsize-$minsize);
		
		if($L>$Lmax){
		    $Lmax = $L;
		    $alphamax = $alpha;
		    $size0max = $size0;
		    $rmax = $r;
		    $varsizemax = $varsize;
		    $vartmax = $vart;
		    $Tendmax = $Tend;
		}
	    }

	    #Get residuals (only if cell lived less than 100 time steps) and cell divided
	    if($Tendmax < 100 && $endtype eq "div"){
		for($t=0;$t<$Tendmax;++$t){
		    $d = $size[$t]-$size0max-$alphamax*($time[$t]-$time[0]);
		    if($avnumfrombottom > 0){
			++$countbegin[$t];
			$resbegin[$t] += $d;
			$tback = $Tendmax-1-$t;
			$endt = 99-$tback;
			++$countend[$endt];
			$resend[$endt] += $d;
		    }
		    else{
			++$countbeginbot[$t];
			$resbeginbot[$t] += $d;
			$tback = $Tendmax-1-$t;
			$endt = 99-$tback;
			++$countendbot[$endt];
			$resendbot[$endt] += $d;
		    }
		}
	    }

	    #Get error bars on the estimated slope and start size
	    if($rmax >= 1.0){
		$rmax = 0.999999;
	    }
	    
	    $sigmaslope = sqrt((1.0-$rmax*$rmax)*$varsizemax/(($Tendmax-1.0)*$vartmax));
	    $sigmasize0 = sqrt($varsizemax*(1.0-$rmax*$rmax)/($Tendmax-1.0));	    
	    #estimate end from the total length (not just fitted length) and growth rate
	    $delt = ($time[$T-1]-$time[0])/($T-1);
	    $size0max -= 0.5*$alpha*$delt;
	    $dsize = $alpha*$delt*$T;
	    $endsize = $musize + $alpha*(0.5*$delt+$time[$T-1]-$time[0]-$mut);
	    $sigmadsize = $sigmaslope*$delt*$T;
	    	    
	    $Tee = $time[$Tendmax-1];
	    #Note, error on edpoint is same as error on start point
	    
	    #Fit fluorescence as a linear function of cell size
	    $startfluo = $fluo[0];
	    $endfluo = $fluo[$T-1];
	    $delfluo = $endfluo-$startfluo;
	    $muf = 0;
	    $musize = 0;
	    $varsize = 0;
	    $varf = 0;
	    $covar = 0;
	    for($t=0;$t<$T;++$t){
		$realsize = exp($size[$t]);
		$musize += $realsize;
		$muf += $fluo[$t];
		$varsize += $realsize*$realsize;
		$varf += $fluo[$t]*$fluo[$t];
		$covar += $realsize*$fluo[$t];
	    }
	    $musize /= $T;
	    $muf /= $T;
	    $covar /= $T;
	    $varsize /= $T;
	    $varf /= $T;
	    $covar  -= $musize*$muf;
	    $varsize -= $musize*$musize;
	    $varf -= $muf*$muf;
	    $rfluo = $covar/sqrt($varsize*$varf);
	    $fluoslope = $covar/$varsize;
	    #Now estimate start and end fluorescence from the fit as a function of cell-size
	    $realstartsize = exp($size0max);
	    $realendsize = exp($endsize);
	    $startfluosize = $muf + $fluoslope*($realstartsize-$musize);
	    $endfluosize = $muf+$fluoslope*($realendsize-$musize);
	    $fluozero = $muf-$fluoslope*$musize;

	    #Take into account that size estimates are noisy?
	    #Put error bar on fluo measurement
	    

	    #Fit fluorescence as a function of time (linear again)
	    $mut = 0;
	    $vart = 0;
	    $covar = 0;
	    for($t=0;$t<$T;++$t){
		$curtime = $time[$t]-$time[0];
		$mut += $curtime;
		$vart += $curtime*$curtime;
		$covar += $fluo[$t]*$curtime;
	    }
	    $mut /= $T;
	    $vart /= $T;
	    $covar /= $T;
	    $vart -= $mut*$mut;
	    $covar -= $mut*$muf;
	    $fluoslopetime = $covar/$vart;
	    $rfluotime = $covar/sqrt($varf*$vart);
	    $startfluotime = $muf-($mut+0.5*$delt)*$fluoslopetime;
	    $endfluotime = $muf + $fluoslopetime*($time[$T-1]-$time[0]+0.5*$delt-$mut);

	    $avnumfrombottom /= $T;
	    $avpos /= $T;

	    	    	    
	    print "$celline $mintop $avpos $avnumfrombottom $Tee $rmax $alphamax $sigmaslope $size0max $endsize $sigmasize0 $dsize $sigmadsize $startfluo $endfluo $rfluo $fluoslope $fluozero $startfluosize $endfluosize $rfluotime $fluoslopetime $startfluotime $endfluotime\n";
	}

	#deal with the new line with cell info
	chomp($line);
	$celline = $line;
	@s = split(/\s+/,$line);
	$endtype = $s[7];
	$celline =~ s/^\>CELL\s+//;
	
	
			
	@time = ();
	@size = ();
	@florsize = ();
	@fluo = ();
	$mintop = 10000;
	#initialize to values that are guaranteed to be beaten
	$maxsize = -1000;
	$minsize = 10000000;
	$avpos = 0;
	$avnumfrombottom = 0;
	$t = 0;
    }
    # frame line (rather than cell line) 
    elsif($_ !~ /^\#/){
	chomp($_);
	@s = split(/\s+/,$_);
	push @time, $s[0];
	++$t;
	$cursize = log($s[5]);
	push @size, $cursize;
	if($cursize>$maxsize){
		$maxsize = $cursize;
	}
	if($cursize<$minsize){
		$minsize = $cursize;
	}
	$curfluo = $s[8];
	push @fluo, $curfluo;

	$curtop = $s[1];
	$curbottom = $s[2];
	$curpos = ($curtop+$curbottom)/2.0;
	$avpos += $curpos;
	if($curtop < $mintop){
	    $mintop = $curtop;
	}

	$curflorsize = log($s[6]);
	push @florsize, $curflorsize;
	$avnumfrombottom += ($s[4]-$s[3]);
    }	
}
# close(G); 

# open(G,$residfile);
print "\n\n\n>RESIDUALS\n";

for($t=0;$t<100;++$t){
    if($countbegin[$t] > 0){
	$resbegin[$t] /= $countbegin[$t];
    }
    if($countend[$t] > 0){
	$resend[$t] /= $countend[$t];
    }
    if($countbeginbot[$t] > 0){
	$resbeginbot[$t] /= $countbeginbot[$t];
    }
    if($countendbot[$t] > 0){
	$resendbot[$t] /= $countendbot[$t];
    }
    
    print $t, "\t", $countbegin[$t], "\t", $resbegin[$t], "\t", $countend[$t], "\t", $resend[$t], "\t", $countbeginbot[$t], "\t", $resbeginbot[$t], "\t", $countendbot[$t], "\t", $resendbot[$t],"\n";
}
# close(G);
