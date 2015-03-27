#!/usr/bin/tclsh

#########################################################
# Tcl script for calculating the INDIVIDUAL free-energy #
# contributions from the restraints, using FEP.         #
#                                                       #
# Usage: ./deltaA_restr.tcl run.dat > A.dat             #
#                                                       #
# Before using, extract the dA/dLambda data from        #
# the log file using:                                   #
# grep "dA/dLambda" run.log > run.dat                   #
#                                                       #
# Author: Yuncheng Mao <catmyc@gmail.com>               #
#              <maoyuncheng@mail.nankai.edu.cn>         #
#                                                       #
#########################################################
# >>>>> Note that the path to tclsh may differ in different machines <<<<<
# >>>>> e.g. tclsh may be located at /usr/local/bin/tclsh            <<<<<

	
if { $argc == 0 } {
	puts "Usage: $argv0 inputfile outputfile"
	exit
} elseif { $argc != 2 } {
	puts "Invalid parameters!"
	puts "Usage: $argv0 inputfile outputfile"
	exit
}

proc deltaA_restr { input output } {
	puts "Opening file $input to read..."
	set infil [open $input r]

	# mark individual restraints
	set ir 1 ; #index of restrants

	set oldL -1

	while { [gets $infil line] != -1 } {
		set newL [lindex $line 2]
		set newGradA [lindex $line 4]
		if { $newL == $oldL } {
			# move to the next restraint
			incr ir
			lappend L($ir) $newL
			lappend gradA($ir) $newGradA
		} else {
			# reset ir
			set ir 1
			lappend L($ir) $newL
			lappend gradA($ir) $newGradA
		}
		set oldL $newL
		#puts "NewL is $newL"
		#puts "NewGradA is $newGradA"
	}
	# Finished reading inputfile.
	puts "$ir restraints recorded."
	puts "Closing file $input..."
	close $infil

	puts "Start calculating the free energy..."
	# Generate a dA, A list
	# first dLambda list
	
	# running through each restraint
	for {set i 1} {$i <= $ir} {incr i} {
		puts "Processing the free energy of Restraint $i"
		foreach L2 [lrange $L($i) 1 end] L1 [lrange $L($i) 0 end-1] {
			lappend dL($i) [expr $L2 - $L1]
		}
		
		# Construct the dA list
		# remove the 0 value of grad_list at lambda=0 first
		set ind [lsearch $L($i) 0]
		set gradA($i) [lreplace $gradA($i) $ind $ind]
		foreach g $gradA($i) dl $dL($i) {
			lappend dA($i) [expr $g * $dl]
		}
		
		# the A list - the free energy profile
		set a 0;
		lappend A($i) $a 
		foreach da $dA($i) {
			set a [expr {$a+$da}]
			lappend A($i) $a
		}

	}

	
	puts "Calculation finished."
	# output 
	puts "Writing to output files..."
	for {set i 1} {$i <= $ir} {incr i} {
		set ofil [open $output.restraint.$i.dat w]
		puts $ofil "# lambda       A         dA"
		set dA($i) [concat 0 $dA($i)]
		foreach l $L($i) a $A($i) da $dA($i) {
			set outstring [format "%-8s  %9.4f  %9.4f" $l $a $da]
			puts $ofil $outstring
		}
		close $ofil
		puts "File $output.restraint.$i.dat written."
	}
	puts "Altogether $ir files created."



}

foreach {inname outname} $argv {}
deltaA_restr $inname $outname

