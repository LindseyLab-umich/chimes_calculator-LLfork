FORCE_FIELDS=../../serial_interface/tests/force_fields
all:
	params_to_jit.py $(FORCE_FIELDS)/published_params.liqC.2b.cubic.txt LIQC2B
	params_to_jit.py $(FORCE_FIELDS)/published_params.liqC.2+3b.cubic.txt LIQC3B
	params_to_jit.py $(FORCE_FIELDS)/published_params.liqCO.2+3b.cubic.txt LIQCO
	params_to_jit.py $(FORCE_FIELDS)/validated_params.CO2400K.2+3+4b.Tersoff.special.offsets.relabel.txt CO_RELABEL
	params_to_jit.py $(FORCE_FIELDS)/published_params.CO2400K.2+3+4b.Tersoff.special.offsets.txt CO2400K
	params_to_jit.py $(FORCE_FIELDS)/published_params.DNTF.2+3b.Tersoff.offsets.DFTB.txt DNTF
	params_to_jit.py $(FORCE_FIELDS)/published_params.HN3.2+3+4b.Tersoff.special.offsets.txt HN3
	params_to_jit.py $(FORCE_FIELDS)/published_params.TATB.txt TATB
	params_to_jit.py $(FORCE_FIELDS)/validated_params.TiO2.2+3b.Tersoff.txt TIO2
	params_to_jit.py $(FORCE_FIELDS)/test_params.CHON.txt CHON
	params_to_jit.py $(FORCE_FIELDS)/test_params.h2o_2bcheby.txt H2O2B
	params_to_jit.py $(FORCE_FIELDS)/test_params.HN3.penalty.txt HN3_PENALTY

clean:
	rm chimesjit.m
