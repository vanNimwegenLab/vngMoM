#!/import/bc2/soft/bin/perl5/perl


$infile = shift(@ARGV);
open(F,$infile) || die "cannot open infile $infile\n";
#get first comment line
$_ = <F>;
my @time = ();
my @size = ();
my @florsize = ();
my @fluo = ();
my $t = 0;


#1:laneID 
#2:cellID  
#3:parID   
#4:daughter-type   
#5:start-time      
#6:end-time        
#7:end-type        
#8:ismintoppixel   
#9:avpos   
#10:avnumfrombottom 
#11:last-time-fit   
#12:Rsize_vs_time   
#13:slope_size_vs_time      
#14:se_slope_size_vs_time   
#15:logsizestart    
#16:logsize_end    
#17:logsize-start-errorbar  
#18:dlogsize        
#19:error_dlogsize  
#20:measfluostart   
#21:measfluoend
#22:Rfluo_vs_size  
#23:slope_fluo_vs_size      
#24:fluo_zero_micron        
#25:est_fluostart_size      
#26:est_fluoend_size        
#27:Rfluo_vs_time   
#28:dslope_fluo_vs_time      
#29:est_fluo_start_time     
#30:est_fluo_end_time       
#31:slope_logfluo_vs_logsize        
#32:se_slope_logfluo_vs_logsize     
#33:slope_fluo_vs_size      
#34:se_slope_fluo_vs_size


if($infile =~ /^([\S\d\/\_]*)parsed\_cellstats\_([\S\_\d]+)\s*$/){
    $prefix = $1;
    $file = $2;
    $outfile = ">" . $prefix . "growth_fluo_" . $file;
}
else{
    $outfile = ">growth_fluo_" . $infile;
}
open(G,$outfile);
print $outfile;
print G "#laneID\tcellID\tparID\tdaughter-type\tstart-time\tend-time\tend-type\tismintoppixel\tavpos\tavnumfrombottom\tlast-time-fit\tRsize_vs_time\tslope_size_vs_time\tse_slope_size_vs_time\tlogsizestart\tlogsize_end\tlogsize-start-errorbar\tdlogsize\terror_dlogsize\tmeasfluostart\tmeasfluoend\tRfluo_vs_size\tslope_fluo_vs_size\tfluo_zero_micron\test_fluostart_size\test_fluoend_size\tRfluo_vs_time\tdsope_fluo_vs_time\test_fluo_start_time\test_fluo_end_time\tslope_logfluo_vs_logsize\tse_slope_logfluo_vs_logsize\tslope_fluo_vs_size\tse_slope_fluo_vs_size\n";

#We only use cells that do not go below pixel $maxheight and that have at least $minlifetime time points in their lifetime!!
$maxheight = 100;
$minlifetime= 4;
while(<F>){
    #new cell
    if($_ =~ /\>CELL/){
	$line = $_;
	#check if there is a previous example to parse
	if($t>=$minlifetime && $mintop >= $maxheight){
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

	    #Get error bars on the estimated slope and start size
	    if($rmax >= 1.0){
		$rmax = 0.999999;
	    }
	    
	    $varslope = (1.0-$rmax*$rmax)*$varsizemax/(($Tendmax-1.0)*$vartmax);
	    if($varslope < 0){
		$varslope = 0.0;
	    }
	    $sigmaslope = sqrt($varslope);
	    $varsize0 = $varsizemax*(1.0-$rmax*$rmax)/($Tendmax-1.0);
	    if($varsize0 < 0){
		$varsize0 = 0;
	    }
	    $sigmasize0 = sqrt($varsize0);
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

	    #Fit log(fluorescence) as a function of log(estimated size from exponential fit against time)
	    #Assume that there is NO OFFSET
	    #Assumes fluorescence is proportional to size + multiplicative noise
	    $muf = 0;
	    $musize = 0;
	    $varsize = 0;
	    $varf = 0;
	    $covar = 0;
	    for($t=0;$t<$T;++$t){
		$estsize = $size0max+$alphamax*($time[$t]-$time[0]);
		$musize += $estsize;
		$logfluo = log($fluo[$t]);
		$muf += $logfluo;
		$varsize += $estsize*$estsize;
		$varf += $logfluo*$logfluo;
		$covar += $estsize*$logfluo;
	    }
	    $musize /= $T;
	    $muf /= $T;
	    $loglam = $muf-$musize;

	    $covar /= $T;
	    $varsize /= $T;
	    $varf /= $T;
	    $covar  -= $musize*$muf;
	    $varsize -= $musize*$musize;
	    $varf -= $muf*$muf;
	    $errloglam = sqrt(($varsize+$varf-2*$covar)/$T);

	    #Fit fluorescence as a function of estimated size from exponential fit against time
	    #Assume that there is NO OFFSET
	    #Assumes fluorescence is proportional to size + constant noise
	    $muf = 0;
	    $musize = 0;
	    $varsize = 0;
	    $varf = 0;
	    $covar = 0;
	    for($t=0;$t<$T;++$t){
		$estsize = exp($size0max+$alphamax*($time[$t]-$time[0]));
		$musize += $estsize;
		$ffll = $fluo[$t];
		$muf += $ffll;
		$varsize += $estsize*$estsize;
		$varf += $ffll*$ffll;
		$covar += $estsize*$ffll;
	    }
	    $musize /= $T;
	    $muf /= $T;
	    $covar /= $T;
	    $varsize /= $T;
	    $varf /= $T;
	    $lam = $covar/$varsize;
	    $errlam = sqrt(($varf/$varsize - $lam*$lam)/$T);
	    
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

	    	    	    
	    print G "$celline $mintop $avpos $avnumfrombottom $Tee $rmax $alphamax $sigmaslope $size0max $endsize $sigmasize0 $dsize $sigmadsize $startfluo $endfluo $rfluo $fluoslope $fluozero $startfluosize $endfluosize $rfluotime $fluoslopetime $startfluotime $endfluotime $loglam $errloglam $lam $errlam\n";
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
close(G); 


