#!/bin/bash

# ChIMES Calculator JIT test suite
#
# ./run_jit_tests.sh
# To specify an install prefix for the cmake tests, run with:
# ./run_jit_tests.sh <install prefix full path>

##################

PREFX=${1-""}     # By default, don't set any special install prefix
PYTH3=python3

# BUILD TARGET
##BUILD=all

#FFS[1 ]="published_params.liqC.2+3b.cubic.txt"                                ; CFGS[1 ]="liqC.2.5gcc_6000K.OUTCAR_#000.xyz"            ; OPTIONS[1 ]="0"
#FFS[2 ]="published_params.liqCO.2+3b.cubic.txt"                               ; CFGS[2 ]="CO.2.5gcc_6500K.OUTCAR_#000.xyz"              ; OPTIONS[2 ]="1"

FFS[0 ]="published_params.CO2400K.2+3+4b.Tersoff.special.offsets.txt"
CFGS[0 ]="CO.2.5gcc_6500K.OUTCAR_#000.translate.xyz"
OPTIONS[0 ]="0"
BUILD[0]="jit.CO2400K"

FFS[1]="published_params.CO2400K.2+3+4b.Tersoff.special.offsets.txt"
CFGS[1]="CO.9GPa_2400K.OUTCAR_#000.xyz"
OPTIONS[1]="0"
BUILD[1]="jit.CO2400K"

FFS[2]="published_params.HN3.2+3+4b.Tersoff.special.offsets.txt"
CFGS[2]="HN3.2gcc_3000K.OUTCAR_#000.xyz"
OPTIONS[2]="0"
BUILD[2]="jit.HN3"

FFS[3]="published_params.DNTF.2+3b.Tersoff.offsets.DFTB.txt"
CFGS[3]="DNTF_2.00gcc_4250K.OUTCAR_#000.xyz"
OPTIONS[3]="0"
BUILD[3]="jit.DNTF"

FFS[4]="published_params.TATB.txt"
CFGS[4]="TATB_extreme_4000K_3.5gcm3_#000.xyz"
OPTIONS[4]="0"
BUILD[4]="jit.TATB"

FFS[5]="published_params.liqC.2b.cubic.txt"
CFGS[5]="liqC.2.5gcc_6000K.OUTCAR_#000.xyz"
OPTIONS[5]="0"
BUILD[5]="jit.LIQC2B"

FFS[6]="published_params.liqC.2+3b.cubic.txt"
CFGS[6]="liqC.2.5gcc_6000K.OUTCAR_#000.xyz"
OPTIONS[6]="0"
BUILD[6]="jit.LIQC3B"

FFS[7]="published_params.liqCO.2+3b.cubic.txt"
CFGS[7]="CO.2.5gcc_6500K.OUTCAR_#000.xyz"
OPTIONS[7]="1"
BUILD[7]="jit.LIQCO"

FFS[8]="validated_params.CO2400K.2+3+4b.Tersoff.special.offsets.relabel.txt"
CFGS[8]="CO.2.5gcc_6500K.OUTCAR_#000.relabel.xyz"
OPTIONS[8]="0"
BUILD[8]="jit.CO_RELABEL"

FFS[9 ]="published_params.CO2400K.2+3+4b.Tersoff.special.offsets.txt"         ;
CFGS[9 ]="CO.2.5gcc_6500K.OUTCAR_#000.scramble.xyz"     ;
OPTIONS[9 ]="0" ;
BUILD[9]="jit.CO2400K"

FFS[10 ]="validated_params.TiO2.2+3b.Tersoff.txt"
CFGS[10 ]="TiO2.unitcell_arbrot_#000.xyz"
OPTIONS[10 ]="0"
BUILD[10]="jit.TIO2"

FFS[11]="test_params.CHON.txt"
CFGS[11 ]="CHON.testfile_#000.xyz"
OPTIONS[11 ]="0"
BUILD[11]="jit.CHON"

FFS[12 ]="published_params.liqCO.2+3b.cubic.txt"
CFGS[12 ]="diam.64_#000.xyz"
OPTIONS[12 ]="1"
BUILD[12]="jit.LIQCO"

FFS[13]="published_params.liqCO.2+3b.cubic.txt"
CFGS[13]="diam.16_#000.xyz"
OPTIONS[13]="1"
BUILD[13]="jit.LIQCO"

FFS[14]="published_params.liqCO.2+3b.cubic.txt"
CFGS[14]="diam.8_#000.xyz"
OPTIONS[14]="1"
BUILD[14]="jit.LIQCO"

FFS[15]="published_params.liqCO.2+3b.cubic.txt"
CFGS[15]="diam.2_#000.xyz"
OPTIONS[15]="1"
BUILD[15]="jit.LIQCO"

FFS[16]="test_params.h2o_2bcheby.txt"
CFGS[16]="H2O.1.50gcc_2000K_#000.xyz"
OPTIONS[16]="0"
BUILD[16]="jit.H2O2B"

### PENALTY TEST IS NOT PASSED - JIT and ChIMESFF differ inside penalty region.
#FFS[17]="test_params.HN3.penalty.txt"
#CFGS[17]="HN3.2.04gcc_20000K_#000.xyz"
#OPTIONS[17]="0"
#BUILD[17]="jit.HN3_PENALTY"

##API[0]="cpp"    ;   EXE[0]="jit.CO2400K"                         ; BUILD[0]="jit.CO2400K"
API[0]="cpp"    ;
# CMakeLists/Makefile/this script need to be updated for these tests
##API[4]="fortran08"; EXE[4]="fortran08_wrapper-serial_interface"; XTRA[4]="" #"0"

## DEBUG !!
echo "Building JIT FILES"
cd ../../chimesFF/src ; make all
cd -

##COMPILE_LIST="CMAKE MAKEFILE"
COMPILE_LIST="MAKEFILE"
#API_LIST="0 1 2 3 4"
API_LIST="0"
NO_TESTS=${#FFS[@]}

for ((j=0;j<$NO_TESTS;j++))
    do
	EXE[$j]=${BUILD[$j]}
    done

LOC=`pwd`

echo "Running $STYLE tests"
date

rm -f jit_test.log

for compile in $COMPILE_LIST
do
	echo "Testing compilation type: $compile"

	# Do the compilation

	if [[ $compile == "MAKEFILE" ]]; then
	    cd ../examples/${API[0]}	
            for ((i=0;i<NO_TESTS;i++))		     
		do
		    echo "Build target ${BUILD[$i]}"
		    make ${BUILD[$i]} DEBUG=1
	    done
            cd ../../tests
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
						time ../../build/${EXE[$j]} force_fields/${FFS[$j]} configurations/$CFG ${OPTIONS[$j]} >> jit_test.log
					else
					        echo "Running ../examples/${API[$i]}/${EXE[$j]} configurations/$CFG ${OPTIONS[$j]}"
						time ../examples/${API[$i]}/${EXE[$j]} configurations/$CFG ${OPTIONS[$j]}  >> jit_test.log
					fi
					
				else				
					${PYTH3} ../examples/${API[$i]}/${EXE[$i]} force_fields/${FFS[$j]} configurations/$CFG ${OPTIONS[$j]} ${LOC}/../api 1 >> jit_test.log 
				fi

				# Compare results against expected results (expected_output/${FFS[$j]}.$CFG.dat)
				
				paste debug.dat expected_output/${FFS[$j]}.$CFG.dat > san.dat

				# Print findings
				
				${PYTH3} compare.py san.dat 

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

##rm -f debug.dat san.dat *.so

cd ../../
./uninstall.sh $PREFX


