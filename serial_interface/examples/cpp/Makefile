
CXX   = g++ -O3 -std=c++11
DEBUG = 0
VERBOSE = 1

CXX_LOC = $(realpath .)
CHIMES_SRC=$(CXX_LOC)/../../../chimesFF/src
SERIAL_SRC=$(CXX_LOC)/../../src

CXX_SRC = main.cpp $(CXX_LOC)/../../src/serial_chimes_interface.cpp $(CHIMES_SRC)/chimesFF.cpp
CXX_HDR = $(CXX_LOC)/../../src/serial_chimes_interface.h $(CHIMES_SRC)/chimesFF.h
JIT_SRC = mainjit.cpp $(CXX_LOC)/../../src/serial_jit_interface.cpp $(CHIMES_SRC)/chimesJIT.cpp \
          $(CHIMES_SRC)/chimesJITev.cpp  $(CHIMES_SRC)/chimesJITcluster.cpp
JIT_CO2400K_SRC = mainjit.cpp $(CXX_LOC)/../../src/serial_jit_interface.cpp $(CHIMES_SRC)/chimesJIT.cpp \
          $(CHIMES_SRC)/chimesJITev.CO2400K.cpp  $(CHIMES_SRC)/chimesJITcluster.CO2400K.cpp
JIT_HN3_SRC = mainjit.cpp $(CXX_LOC)/../../src/serial_jit_interface.cpp $(CHIMES_SRC)/chimesJIT.cpp \
          $(CHIMES_SRC)/chimesJITev.HN3.cpp  $(CHIMES_SRC)/chimesJITcluster.HN3.cpp
JIT_DNTF_SRC = mainjit.cpp $(CXX_LOC)/../../src/serial_jit_interface.cpp $(CHIMES_SRC)/chimesJIT.cpp \
          $(CHIMES_SRC)/chimesJITev.DNTF.cpp  $(CHIMES_SRC)/chimesJITcluster.DNTF.cpp
JIT_TATB_SRC = mainjit.cpp $(CXX_LOC)/../../src/serial_jit_interface.cpp $(CHIMES_SRC)/chimesJIT.cpp \
          $(CHIMES_SRC)/chimesJITev.TATB.cpp  $(CHIMES_SRC)/chimesJITcluster.TATB.cpp
JIT_LIQC2B_SRC = mainjit.cpp $(CXX_LOC)/../../src/serial_jit_interface.cpp $(CHIMES_SRC)/chimesJIT.cpp \
          $(CHIMES_SRC)/chimesJITev.LIQC2B.cpp  $(CHIMES_SRC)/chimesJITcluster.LIQC2B.cpp
JIT_LIQC3B_SRC = mainjit.cpp $(CXX_LOC)/../../src/serial_jit_interface.cpp $(CHIMES_SRC)/chimesJIT.cpp \
          $(CHIMES_SRC)/chimesJITev.LIQC3B.cpp  $(CHIMES_SRC)/chimesJITcluster.LIQC3B.cpp
JIT_LIQCO_SRC = mainjit.cpp $(CXX_LOC)/../../src/serial_jit_interface.cpp $(CHIMES_SRC)/chimesJIT.cpp \
          $(CHIMES_SRC)/chimesJITev.LIQCO.cpp  $(CHIMES_SRC)/chimesJITcluster.LIQCO.cpp
JIT_CO_RELABEL_SRC = mainjit.cpp $(CXX_LOC)/../../src/serial_jit_interface.cpp $(CHIMES_SRC)/chimesJIT.cpp \
          $(CHIMES_SRC)/chimesJITev.CO_RELABEL.cpp  $(CHIMES_SRC)/chimesJITcluster.CO_RELABEL.cpp
JIT_TIO2_SRC = mainjit.cpp $(CXX_LOC)/../../src/serial_jit_interface.cpp $(CHIMES_SRC)/chimesJIT.cpp \
          $(CHIMES_SRC)/chimesJITev.TIO2.cpp  $(CHIMES_SRC)/chimesJITcluster.TIO2.cpp
JIT_CHON_SRC = mainjit.cpp $(CXX_LOC)/../../src/serial_jit_interface.cpp $(CHIMES_SRC)/chimesJIT.cpp \
          $(CHIMES_SRC)/chimesJITev.CHON.cpp  $(CHIMES_SRC)/chimesJITcluster.CHON.cpp
JIT_H2O2B_SRC = mainjit.cpp $(CXX_LOC)/../../src/serial_jit_interface.cpp $(CHIMES_SRC)/chimesJIT.cpp \
          $(CHIMES_SRC)/chimesJITev.H2O2B.cpp  $(CHIMES_SRC)/chimesJITcluster.H2O2B.cpp
JIT_HN3_PENALTY_SRC = mainjit.cpp $(CXX_LOC)/../../src/serial_jit_interface.cpp $(CHIMES_SRC)/chimesJIT.cpp \
          $(CHIMES_SRC)/chimesJITev.HN3_PENALTY.cpp  $(CHIMES_SRC)/chimesJITcluster.HN3_PENALTY.cpp
##
JIT_HDR = $(CXX_LOC)/../../src/serial_jit_interface.h $(CHIMES_SRC)/chimesJIT.h
JIT_HN3_HDR=$(CHIMES_SRC)/chimesJITev.HN3.h
JIT_CO2400K_HDR=$(CHIMES_SRC)/chimesJITev.CO2400K.h
JIT_DNTF_HDR=$(CHIMES_SRC)/chimesJITev.DNTF.h
JIT_TATB_HDR=$(CHIMES_SRC)/chimesJITev.TATB.h
JIT_LIQC2B_HDR=$(CHIMES_SRC)/chimesJITev.LIQC2B.h
JIT_LIQC3B_HDR=$(CHIMES_SRC)/chimesJITev.LIQC3B.h
JIT_LIQCO_HDR=$(CHIMES_SRC)/chimesJITev.LIQCO.h
JIT_CO_RELABEL_HDR=$(CHIMES_SRC)/chimesJITev.CO_RELABEL.h
JIT_TIO2_HDR=$(CHIMES_SRC)/chimesJITev.TIO2.h
JIT_CHON_HDR=$(CHIMES_SRC)/chimesJITev.CHON.h
JIT_H2O2B_HDR=$(CHIMES_SRC)/chimesJITev.H2O2B.h
JIT_HN3_PENALTY_HDR=$(CHIMES_SRC)/chimesJITev.HN3_PENALTY.h

clean:

clean-all:
	rm -f CPP-interface jit.*

all: $(CXX_SRC) $(CXX_HDR)
	$(CXX) $(CXX_SRC)  -o CPP-interface -I $(SERIAL_SRC) -I $(CHIMES_SRC)/ -DDEBUG=${DEBUG} -DVERBOSE=${VERBOSE}

jit.CO2400K : $(JIT_CO2400K_SRC) $(JIT_HDR) $(JIT_CO2400K_HDR)
	cp $(CHIMES_SRC)/chimesJITev.CO2400K.h $(CHIMES_SRC)/chimesJITev.h
	$(CXX) $(JIT_CO2400K_SRC)  -o jit.CO2400K -I $(SERIAL_SRC) -I $(CHIMES_SRC)/ -DDEBUG=${DEBUG} -DVERBOSE=${VERBOSE}

jit.HN3 : $(JIT_HN3_SRC) $(JIT_HDR) $(JIT_HN3_HDR)
	cp $(CHIMES_SRC)/chimesJITev.HN3.h $(CHIMES_SRC)/chimesJITev.h
	$(CXX) $(JIT_HN3_SRC)  -o jit.HN3 -I $(SERIAL_SRC) -I $(CHIMES_SRC)/ -DDEBUG=${DEBUG} -DVERBOSE=${VERBOSE}

jit.DNTF : $(JIT_DNTF_SRC) $(JIT_HDR) $(JIT_DNTF_HDR)
	cp $(CHIMES_SRC)/chimesJITev.DNTF.h $(CHIMES_SRC)/chimesJITev.h
	$(CXX) $(JIT_DNTF_SRC)  -o jit.DNTF -I $(SERIAL_SRC) -I $(CHIMES_SRC)/ -DDEBUG=${DEBUG} -DVERBOSE=${VERBOSE}

jit.TATB : $(JIT_TATB_SRC) $(JIT_HDR) $(JIT_TATB_HDR)
	cp $(CHIMES_SRC)/chimesJITev.TATB.h $(CHIMES_SRC)/chimesJITev.h
	$(CXX) $(JIT_TATB_SRC)  -o jit.TATB -I $(SERIAL_SRC) -I $(CHIMES_SRC)/ -DDEBUG=${DEBUG} -DVERBOSE=${VERBOSE}

jit.LIQC2B : $(JIT_LIQC2B_SRC) $(JIT_HDR) $(JIT_LIQC2B_HDR)
	cp $(CHIMES_SRC)/chimesJITev.LIQC2B.h $(CHIMES_SRC)/chimesJITev.h
	$(CXX) $(JIT_LIQC2B_SRC)  -o jit.LIQC2B -I $(SERIAL_SRC) -I $(CHIMES_SRC)/ -DDEBUG=${DEBUG} -DVERBOSE=${VERBOSE}

jit.LIQC3B : $(JIT_LIQC3B_SRC) $(JIT_HDR) $(JIT_LIQC3B_HDR)
	cp $(CHIMES_SRC)/chimesJITev.LIQC3B.h $(CHIMES_SRC)/chimesJITev.h
	$(CXX) $(JIT_LIQC3B_SRC)  -o jit.LIQC3B -I $(SERIAL_SRC) -I $(CHIMES_SRC)/ -DDEBUG=${DEBUG} -DVERBOSE=${VERBOSE}

jit.LIQCO : $(JIT_LIQCO_SRC) $(JIT_HDR) $(JIT_LIQCO_HDR)
	cp $(CHIMES_SRC)/chimesJITev.LIQCO.h $(CHIMES_SRC)/chimesJITev.h
	$(CXX) $(JIT_LIQCO_SRC)  -o jit.LIQCO -I $(SERIAL_SRC) -I $(CHIMES_SRC)/ -DDEBUG=${DEBUG} -DVERBOSE=${VERBOSE}

jit.CO_RELABEL : $(JIT_CO_RELABEL_SRC) $(JIT_HDR) $(JIT_CO_RELABEL_HDR)
	cp $(CHIMES_SRC)/chimesJITev.CO_RELABEL.h $(CHIMES_SRC)/chimesJITev.h
	$(CXX) $(JIT_CO_RELABEL_SRC)  -o jit.CO_RELABEL -I $(SERIAL_SRC) -I $(CHIMES_SRC)/ -DDEBUG=${DEBUG} -DVERBOSE=${VERBOSE}

jit.TIO2 : $(JIT_TIO2_SRC) $(JIT_HDR) $(JIT_TIO2_HDR)
	cp $(CHIMES_SRC)/chimesJITev.TIO2.h $(CHIMES_SRC)/chimesJITev.h
	$(CXX) $(JIT_TIO2_SRC)  -o jit.TIO2 -I $(SERIAL_SRC) -I $(CHIMES_SRC)/ -DDEBUG=${DEBUG} -DVERBOSE=${VERBOSE}

jit.CHON : $(JIT_CHON_SRC) $(JIT_HDR) $(JIT_CHON_HDR)
	cp $(CHIMES_SRC)/chimesJITev.CHON.h $(CHIMES_SRC)/chimesJITev.h
	$(CXX) $(JIT_CHON_SRC)  -o jit.CHON -I $(SERIAL_SRC) -I $(CHIMES_SRC)/ -DDEBUG=${DEBUG} -DVERBOSE=${VERBOSE}

jit.H2O2B : $(JIT_H2O2B_SRC) $(JIT_HDR) $(JIT_H2O2B_HDR)
	cp $(CHIMES_SRC)/chimesJITev.H2O2B.h $(CHIMES_SRC)/chimesJITev.h
	$(CXX) $(JIT_H2O2B_SRC)  -o jit.H2O2B -I $(SERIAL_SRC) -I $(CHIMES_SRC)/ -DDEBUG=${DEBUG} -DVERBOSE=${VERBOSE}

jit.HN3_PENALTY : $(JIT_HN3_PENALTY_SRC) $(JIT_HDR) $(JIT_HN3_PENALTY_HDR)
	cp $(CHIMES_SRC)/chimesJITev.HN3_PENALTY.h $(CHIMES_SRC)/chimesJITev.h
	$(CXX) $(JIT_HN3_PENALTY_SRC)  -o jit.HN3_PENALTY -I $(SERIAL_SRC) -I $(CHIMES_SRC)/ -DDEBUG=${DEBUG} -DVERBOSE=${VERBOSE}


