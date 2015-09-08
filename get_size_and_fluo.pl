#!/import/bc2/soft/bin/perl5/perl


### PARAMETERS TO TRANSFORM MEASUREMENTS TO PHYSICAL QUANTITIES ####

#Number of seconds per time point
$delt = 180;
#Micrometers per pixel
$micronperpixel = 0.065;
#GFP per fluorescence intensity unit
$GFPperfluo = 1.0; # 1.0/4.02;


$numargs = @ARGV;
if($numargs != 2){
    die "USAGE new_get_size_and_fluo.pl infile debug\nWhere debug is a binary variable.";
}

$infile = shift(@ARGV);

#If debug is set to one we will print out extra info: rowavs, fluo distributions, and cell_edges
$debug = shift(@ARGV);

if($infile =~ /pos(\d+)/){
		$pos = $1;
}
else{
		die "cannot get position from $infile\n";
}
if($infile =~ /GL(\d+)/){
		$GL = $1;
}
else{
		die "cannot get growth lane id GL from $infile.\n";
}
$laneID = "pos_" . $pos . "_GL_" . $GL;

if($debug){
		if($infile =~ /ExportedCellStats\_(\d+)\_/){
				$date = $1;
		}
		else{
				die "cannot get date from $infile.\n";
		}

		if($infile =~ /^([\S\/\d\_]*)ExportedCellStats/){
				$prefix = $1;
		}
		else{
				die "cannot get prefix from $infile\n";
		}

    $outfile = ">" . $prefix . "cell_edges_" . $date . "_pos" . $pos . "_GL" . $GL;
    open(H,$outfile);
}



#Get row sums for all time points of a cell.
#For each time point:
#Calculate also derivatives of the rowsums (i.e. rowsums[i+1]-rowsums[i])
#If cell is not the bottom cell:
#Get 10 consecutive points with lowest average sum of rowsums. This is assumed to be in the middle of the cell.
#Go left and right from this middle point and find lowest derivative to the left and highest derivative to the right.
#Size is distance between right and left point.
#If bottom cell:
#Find 10 consecutive points with lowest sum of derivative^2, i.e. flattest region
#Consider this in the middle of the cell.
#Go left to find region with lowest derivative
#If rowsum at end is lower than rowsum at end-3: assume cell is pushed into the bottom -> cell end is most negative derivative
#otherwise look for most positive derivative on right.


open(F,$infile) || die "cannot open $infile \n";
print "#CELL lane_ID\tcell_ID\tparent_ID\tdaughter-type\tfirst-time(sec)\tlast-time(sec)\ttype_of_end\n";
print "#time(sec)\tvertical_top\tvertical_bottom\tcell_num_in_lane\ttotal_cell_in_lane\theight(micrometer)\tflorian_height\tfluo_background\tfluo_amplitude\tfluo_middle\tfluo_width\n";
$curid = -1;
my @curstats = ();
my @rs = ();

while(<F>){
    #get line
    $line = $_;
    #remove MAC or WINDOWS carriage returns and replace with newlines
    $line =~ s/\r+//g;
    $line =~ s/\n//;
    #Check if new cell ID
    if($line =~ /id\=(\d+)\;/){
	#Get new ID, parent ID, birth frame, and daughter type
	$newid = $1;
		
	if($line =~ /birth\_frame\=([\-\d]+)/){
            $newbf = $1;
        }
        else{
            die "cannot determine birth_frame from $line\n";
        }
	if($line =~ /pid\=([\-\d]+)/){
	    $newpid = $1;
	}
	else{
	    die "cannot determine parent ID from $line\n";
	}
	if($line =~ /daughter\_type\=(\S+)\s*$/){
	    $newdaughtertype = $1;
	}
	else{
	    die "cannot get daughter type from $line\n";
	}
	
	#print out all information of the previous cell
	if($curid >= 0){
	    $birthtime = $birthframe * $delt;
	    $curtime = $curframe * $delt;
	    print ">CELL $laneID\t$curid\t$pid\t$daughtertype\t$birthtime\t$curtime\t$curend\n";
	    $time = @curstats;
	    for($t=0;$t<$time;++$t){
		print $curstats[$t], "\n";
	    }

	    #print rowsums if debugging
	    if($debug){
		$curoutfile = ">" . $prefix . "rowavs_" . $laneID . "_cell_" . $curid . "_" . $curend . "_t_" . $time;
		open(J,$curoutfile);
		for($row=0;$row<$maxsize;++$row){
		    print J $row;
		    for($tt=0;$tt<$time;++$tt){
			if(defined($rs[$tt][$row])){
			    print J "\t" . $rs[$tt][$row];
			}
			else{
			    print J "\t0";
			}
		    }
		    print J "\n";
		}
		close(J);
	    }
	}
	
	#Now assign the new ID, par ID etc to the current cell.
	$curid = $newid;
	$birthframe = $newbf;
	if($birthframe < 0){
	    $birthframe = 0;
	}
	$pid = $newpid;
	$daughtertype = $newdaughtertype;
	
	#print STDERR "doing $curid with parent $pid birthframe $newbf\n";	
	@curstats = ();
	if($debug){
	    @rs = ();
	    $t = 0;
	    $maxsize = 0;
	}
    }
    #Line with all contrast pixel intensities. Use them to get row averages, and then estimate current size
    elsif($line =~ /ch\=0\;\s+output\=PIXEL\_INTENSITIES\;([\d\.\;\s]+)/){
	$vals = $1;
	$vals =~ s/\;/ /g;
	my @in = split(/\s+/,$vals);
	#Get the row sums
	$size = @in/20;
	if($debug){
	    if($size > $maxsize){
		$maxsize = $size;
	    }
	}
	my @rowsums = ();
	for($row=0;$row<$size;++$row){
	    $rowsum = 0;
	    for($col=0;$col<20;++$col){
		$rowsum += $in[20*$row+$col];
	    }
	    $rowsums[$row] = $rowsum;
	    if($debug){
		$rs[$t][$row] = $rowsum;
	    }
	}
	if($debug){
	    ++$t;
	}

	#Get the derivative of the rowsums
	my @curderiv = ();
	#get derivatives
	for($i=0;$i<$size-1;++$i){
	    $der = $rowsums[$i+1]-$rowsums[$i];
	    push @curderiv, $der;
	}

	#Check if it is the bottom cell
	$isbottom = 0;
	if($cellnum == $cellinlane){
	    $isbottom = 1;
	}
	
	if($isbottom){
	    #find flattest region
	    $minvar = 100000000.0;
	    for($i=0;$i<$size-11;++$i){
		$curvar = 0;
		for($j=0;$j<10;++$j){
		    $curvar += $curderiv[$i+$j]*$curderiv[$i+$j];
		}
		if($curvar < $minvar){
		    $minvar = $curvar;
		    $minpos = $i;
		}
	    }

	    #Get average absolute value in this flattest region (not used right now)
	    $absmid = 0;
	    for($j=0;$j<10;++$j){
		$absmid += $rowsums[$minpos+$j];
	    }
	    $absmid /= 10;
	    
	    $mid = $minpos + 5;

	    #Left end is where derivative is most negative
	    $minderiv = 10000000.0;
	    for($i=0;$i<$mid;++$i){
		if($curderiv[$i] < $minderiv){
		    $minderiv = $curderiv[$i];
		    $imin = $i;
		}
	    }

	    $crit1 = 0;
	    $crit2 = 0;
	    #First criterium: Is the endpoint above the average value in the middle ?
	    if($rowsums[$size-1]  > $absmid){
		$crit1 = 1;
	    }
	    #Second criterium: Is the last data point below the value 3 points back? (ends with decrease)
	    #This is the one we are using right now.
	    if($rowsumes[$size-1] > $rowsums[$size-4]){
		$crit2 = 1;
	    }
	    	    
	    #Right end side. See if there is a final increase (go back 3 points)
	    if($crit2){
		#End went up: Look for the maximum positive derivative (like normal cells)
		$maxderiv = -10000000.0;
		for($i=$mid;$i<$size-1;++$i){
		    if($curderiv[$i] > $maxderiv){
			$maxderiv = $curderiv[$i];
			$imax = $i;
		    }
		}
	    }#Else we have cell pushed to the bottom and look for most negative derivative. In the last 5 points only!
	    else{
		$minderiv = 10000000.0;
		for($i=$size-6;$i<$size-1;++$i){
		    if($curderiv[$i] < $minderiv){
			$minderiv = $curderiv[$i];
			$imax = $i;
		    }
		}
	    }
	}#Else we have a non-bottom cell
	else{
	    #Middle is where 10 consecutive rows have lowest sum of rowsums
	    $minsum = 100000000.0;
	    for($i=0;$i<$size-10;++$i){
		$cursum = 0;
		for($j=0;$j<10;++$j){
		    $cursum += $rowsums[$i+$j];
		}
		if($cursum < $minsum){
		    $minsum = $cursum;
		    $minpos = $i;
		}
	    }
	    $mid = $minpos+5;

	    #Left end is where derivative is most negative
	    $minderiv = 10000000.0;
	    for($i=0;$i<$mid;++$i){
		if($curderiv[$i] < $minderiv){
		    $minderiv = $curderiv[$i];
		    $imin = $i;
		}
	    }

	    #Right end is most positive  derivative
	    $maxderiv = -10000000.0;
	    for($i=$mid;$i<$size-1;++$i){
		if($curderiv[$i] > $maxderiv){
		    $maxderiv = $curderiv[$i];
		    $imax = $i;
		}
	    }
	}
	    
	
	#Estimate of the size
	$mysize = $imax-$imin;
	
	if($debug){
	    print H "ID $curid frame $curframe mysize $mysize mid $mid imax $imax imin $imin maxderiv $maxderiv minderiv $minderiv bottom $isbottom $crit1 $crit2\n"
	}
    }
    #frame line with cell height and position in growth lane.
    elsif($line =~ /frame\=(\d+)/){
	$curframe = $1;
	if($line =~ /height\=([\d\.]+)/){
	    $florheight = $1;
	}
	else{
	    die "cannot get height from: $line\n";
	}

	if($line =~ /pixel\_limits\=\[(\d+)\,(\d+)\]/){
	    $toppixel = $1;
	    $botpixel = $2;
	}
	else{
	    die "cannot get pixel limits from $line\n";
	}

	if($line =~ /pos\_in\_GL\=\[(\d+)\,(\d+)\]/){
	    $cellnum = $1;
	    $cellinlane = $2;
	}
	else{
	    die "cannot get cell number in lane from $line\n";
	}
    }
    #line with fluorescence column intensities
    elsif($line =~ /ch\=1\;\s+output\=COLUMN\_INTENSITIES\;\s+([\d\.\;\s]+)/){
	$vals = $1;
	$vals =~ s/\;/ /g;
	my @in = split(/\s+/,$vals);
	$num = @in;

	#get maximum and minimum
	$max = -10000000;
	$min = 100000000;
	for($k=0;$k<$num;++$k){
	    if($in[$k] > $max){
		$max = $in[$k];
	    }
	    if($in[$k] < $min){
		$min = $in[$k];
	    }
	}
	
	#Initial guess at parameters
	$mu = 50;
	$del= 5.5;
	$diff = 1;
	$b = $min;
	$A = $max-$min;
	my @rho = ();
	while($diff > 0.01){
	    #Set rho, A, and B using EM equation
	    $rhosum = 0;
	    for($i=0;$i<$num;++$i){
		#calc 1/(1+(x-mu)^2/del^2)
		$rho[$i] = ($i-$mu)/$del;
		$rho[$i] = $rho[$i]*$rho[$i];
		$rho[$i] = 1.0/(1.0+$rho[$i]);
		$rhosum += $rho[$i];
	    }
    
	    $bnew = 0;
	    for($i=0;$i<$num;++$i){
		$bnew += $in[$i]*$b/($b+$A*$rho[$i]);
	    }
	    $bnew /= $num;

	    $Anew = 0;
	    for($i=0;$i<$num;++$i){
		$Anew += $in[$i]*$A*$rho[$i]/($b+$A*$rho[$i]);
	    }
	    $Anew /= $rhosum;

	    #Get the new mu using binary search
	    $curmin = 40;
	    $curmax = 60;
	    while($curmax-$curmin > 0.01){
		$curmu = ($curmin+$curmax)/2;
		$deriv = 0;
		for($i=0;$i<$num;++$i){
		    $thisrho = ($i-$curmu)/$del;
		    $thisrho = $thisrho*$thisrho;
		    $thisrho = 1.0/(1.0+$thisrho);
		    $deriv += ($i-$mu)*$thisrho*$thisrho*(-1.0+$in[$i]/($bnew+$Anew*$thisrho));
		}
		if($deriv>0){
		    $curmin = $curmu;
		}
		else{
		    $curmax = $curmu;
		}
	    }
	    $munew = ($curmax+$curmin)/2;

	    #Get new delta (width of Cauchy peak) using binary search
	    $deltamax = 12.0;
	    $deltamin = 3.0;
	    while($deltamax-$deltamin > 0.01){
		$curdel = ($deltamax+$deltamin)/2.0;
		$deriv = 0;
		for($i=0;$i<$num;++$i){
		    $thisrho = ($i-$munew)/$curdel;
		    $thisrho = $thisrho*$thisrho;
		    $thisrho = 1.0/(1.0+$thisrho);
		    $cur = $thisrho*(1.0-$thisrho)*(-1.0 + $in[$i]/($bnew+$Anew*$thisrho));
		    $deriv += $cur;
		}
		if($deriv > 0){
		    $deltamin = $curdel;
		}
		else{
		    $deltamax = $curdel;
		}
	    }
	    $deltanew = $curdel;
	     
	    $diff = 0;
	    $diff += abs($Anew-$A)/($Anew+$A);
	    $diff += abs($bnew-$b)/($b+$bnew);
	    $diff += abs($munew-$mu)/($munew+$mu);
	    $diff += abs($deltanew-$del)/($deltanew+$del);
	    
	    $b = $bnew;
	    $A = $Anew;
	    $mu = $munew;
	    $del = $deltanew;
	}

	if($debug){
	    #print observed and fitted fluorescence intensities
	    $curh = $toppixel-$botpixel+1;
	    $curbox = 20*$curh;
	    $file = ">fluo_id" . $curid . "_frame_" . $curframe . "_curh" . $curh . "_box" . $curbox;
	    open(FF,$file);
	    for($j=0;$j<$num;++$j){
		$thisrho = ($j-$mu)/$del;
		$thisrho = $thisrho*$thisrho;
		$thisrho = 1.0/(1.0+$thisrho);
		$theo = $b + $A*$thisrho;
		print FF $in[$j], " ", $theo, "\n";
	    }
	    close(FF);
	}


	#Since these are last stats for this frame, record stats on this frame.	
	
	#Stats we are going to print:
	#Current frame
	#top and bottom pixels of box
	#Cell number of channel
	#Total cells in channel
	#My height
	#Florian's height
	#Fluoresence data:
	#background
	#amplitude
	#middle (of 100)
	#width
	$curtime = $curframe * $delt;
	$realsize = $mysize * $micronperpixel;
	$florsize = $florheight * $micronperpixel;
	$bGFP = $b * $GFPperfluo;
	$aGFP = $A * $GFPperfluo;
		
	$stats = $curtime . "\t" . $toppixel . "\t" . $botpixel . "\t" . $cellnum . "\t" . $cellinlane . "\t" . $realsize . "\t" . $florsize . "\t" . $bGFP . "\t" . $aGFP . "\t" . $mu . "\t" . $del;
	push @curstats, $stats;
    }
    elsif($line =~ /EXIT/){
	$curend = "exit";
    }
    elsif($line =~ /DIVISION/){
	$curend = "div";
    }
    elsif($line =~ /ENDOFDATA/){
	$curend = "eod";
    }
}
close(F);


# close(G);
close(H);
