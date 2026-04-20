#!/bin/bash

# ChIMES Calculator test suite
#
# For (relatively) quick and dirty testing, run with:
# ./run_tests.sh SHORT
# Otherwise, run with:
# ./run_tests.sh
# To specify an install prefix for the cmake tests, run with:
# ./run_tests.sh <SHORT or LONG> <install prefix full path>


##################

STYLE=${1-"LONG"} # By default,run "LONG" test, but if user runs with "./run_tests SHORT, runs short tests
PREFX=${2-""}     # By default, don't set any special install prefix
PYTH3=python3


FFS[0 ]="published_params.CO2400K.2+3+4b.Tersoff.special.offsets.txt"         ; CFGS[0 ]="CO.2.5gcc_6500K.OUTCAR_#000.translate.xyz"    ; OPTIONS[0 ]="0"
FFS[1 ]="validated_params.CO2400K.2+3+4b.Tersoff.special.offsets.relabel.txt" ; CFGS[1 ]="CO.2.5gcc_6500K.OUTCAR_#000.relabel.xyz"      ; OPTIONS[1 ]="0"
FFS[2]="published_params.CO2400K.2+3+4b.Tersoff.special.offsets.txt"          ; CFGS[2]="CO.9GPa_2400K.OUTCAR_#000.xyz"                 ; OPTIONS[2 ]="0"
FFS[3 ]="published_params.HN3.2+3+4b.Tersoff.special.offsets.txt"             ; CFGS[3 ]="HN3.2gcc_3000K.OUTCAR_#000.xyz"               ; OPTIONS[3 ]="0"
FFS[4 ]="published_params.DNTF.2+3b.Tersoff.offsets.DFTB.txt"                 ; CFGS[4 ]="DNTF_2.00gcc_4250K.OUTCAR_#000.xyz"           ; OPTIONS[4 ]="0"
FFS[5 ]="published_params.TATB.txt"                                           ; CFGS[5 ]="TATB_extreme_4000K_3.5gcm3_#000.xyz"          ; OPTIONS[5 ]="0"
FFS[6 ]="published_params.liqC.2b.cubic.txt"                                  ; CFGS[6 ]="liqC.2.5gcc_6000K.OUTCAR_#000.xyz"            ; OPTIONS[6 ]="0"
FFS[7 ]="published_params.liqC.2+3b.cubic.txt"                                ; CFGS[7 ]="liqC.2.5gcc_6000K.OUTCAR_#000.xyz"            ; OPTIONS[7 ]="0"
FFS[8 ]="published_params.liqCO.2+3b.cubic.txt"                               ; CFGS[8 ]="CO.2.5gcc_6500K.OUTCAR_#000.xyz"              ; OPTIONS[8 ]="1" 
FFS[9 ]="published_params.CO2400K.2+3+4b.Tersoff.special.offsets.txt"         ; CFGS[9 ]="CO.2.5gcc_6500K.OUTCAR_#000.scramble.xyz"     ; OPTIONS[9 ]="0"
FFS[10]="validated_params.TiO2.2+3b.Tersoff.txt"                              ; CFGS[10 ]="TiO2.unitcell_arbrot_#000.xyz"               ; OPTIONS[10]="0"
FFS[11]="test_params.CHON.txt"                                                ; CFGS[11 ]="CHON.testfile_#000.xyz"                      ; OPTIONS[11]="0"
FFS[12]="published_params.liqCO.2+3b.cubic.txt"                               ; CFGS[12 ]="diam.64_#000.xyz"                            ; OPTIONS[12]="1"
FFS[13]="published_params.liqCO.2+3b.cubic.txt"                               ; CFGS[13]="diam.16_#000.xyz"                             ; OPTIONS[13]="1"
FFS[14]="published_params.liqCO.2+3b.cubic.txt"                               ; CFGS[14]="diam.8_#000.xyz"                              ; OPTIONS[14]="1"
FFS[15]="published_params.liqCO.2+3b.cubic.txt"                               ; CFGS[15]="diam.2_#000.xyz"                              ; OPTIONS[15]="1"
FFS[16]="test_params.h2o_2bcheby.txt"                                         ; CFGS[16]="H2O.1.50gcc_2000K_#000.xyz"                   ; OPTIONS[16]="0"
FFS[17]="test_params.HN3.penalty.txt"                                         ; CFGS[17]="HN3.2.04gcc_20000K_#000.xyz"                  ; OPTIONS[17]="0"

API[0]="cpp"    ;   EXE[0]="CPP-interface"                     ; XTRA[0]="" #"2"
API[1]="c"      ;   EXE[1]="C_wrapper-serial_interface"        ; XTRA[1]="" #"2"
API[2]="fortran";   EXE[2]="fortran_wrapper-serial_interface"  ; XTRA[2]="" #"2"
API[3]="python" ;   EXE[3]="main.py"                           ; XTRA[3]="" #"2 1"
# CMakeLists/Makefile/this script need to be updated for these tests
##API[4]="fortran08"; EXE[4]="fortran08_wrapper-serial_interface"; XTRA[4]="" #"0"


COMPILE_LIST="CMAKE MAKEFILE"
##COMPILE_LIST="MAKEFILE"
API_LIST="0 1 2 3"
NO_TESTS=${#FFS[@]}
LOC=`pwd`

echo "Running $STYLE tests"
date

rm -f test.log

for compile in $COMPILE_LIST
do
	echo "Testing compilation type: $compile"

	# Do the compilation

	if [[ $compile == "MAKEFILE" ]]; then
	
		for i in $API_LIST # Cycle through APIs
		do
			cd ../examples/${API[$i]}	
	
			echo "Compiling for API ${API[$i]}"
			echo ""

			if [[ "${API[$i]}" != "python" ]] ; then
	
				make all DEBUG=1
			else
				make all
			fi
			
			cd ../../tests
			
		done		
	
	elif [[ $compile == "CMAKE" ]] ; then
		
		cd ../../
		./install.sh 0 "$PREFX" 1 1 # Set the verbose and LLNL flags true
		cp build/lib-C_wrapper-serial_interface.so  serial_interface/examples/python		
		cd -  

	else
		echo "Error: Unknown compilation method $compile"
		echo "Acceptable values are MAKEFILE and CMAKE"
		echo "Check logic in run_test.sh"
		exit 0
	fi
	
	# Run the tasks
		
	for i in $API_LIST # Cycle through APIs
	do	
		
		idx=1

		for ((j=0;j<NO_TESTS;j++))
		do
			echo "Working on Test $idx of $NO_TESTS for API ${API[$i]}"
			
			for ((k=0; k<10; k++))
			do
				CFG=${CFGS[$j]}
				CFG_PREFIX="${CFG%%_#*}_#"
				CFG_SUFFIX=${CFG##*000}
				
				CFG=${CFG_PREFIX}00${k}${CFG_SUFFIX}

				if [ ! -f configurations/$CFG ] ; then
				    continue
				fi
				
				echo "		...Running $CFG"
				
				# Run the test
				
				if [[ "${API[$i]}" != "python" ]] ; then
				
					if [[ $compile == "CMAKE" ]] ; then
					    time ../../build/${EXE[$i]} force_fields/${FFS[$j]} configurations/$CFG ${OPTIONS[$j]} >> test.log
					else
					    time ../examples/${API[$i]}/${EXE[$i]} force_fields/${FFS[$j]} configurations/$CFG ${OPTIONS[$j]}  >> test.log
					fi
					
				else				
					${PYTH3} ../examples/${API[$i]}/${EXE[$i]} force_fields/${FFS[$j]} configurations/$CFG ${OPTIONS[$j]} ${LOC}/../api 1 >> test.log 2>> test.log 
				fi

				# Compare results against expected results (expected_output/${FFS[$j]}.$CFG.dat)
				
				paste debug.dat expected_output/${FFS[$j]}.$CFG.dat > san.dat

				# Print findings
				
				${PYTH3} compare.py san.dat 

				if [[ $STYLE == "SHORT" ]] ; then
					break
				fi
				
			done
			
			echo "	Test $idx of $NO_TESTS for API ${API[$i]} complete"

			let idx=idx+1
			
		done

	done

done

# Clean up

for i in $API_LIST # Cycle through APIs
do
	cd ../examples/${API[$i]}
	make clean
	rm -f *.so *.a
	cd ../../tests
done

rm -f debug.dat san.dat *.so

cd ../../
./uninstall.sh $PREFX


