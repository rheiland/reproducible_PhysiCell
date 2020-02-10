VERSION := $(shell grep . VERSION.txt | cut -f1 -d:)
PROGRAM_NAME := test_cell_div_bdy

CC := g++
# CC := g++-mp-7 # typical macports compiler name
# CC := g++-7 # typical homebrew compiler name 

# Check for environment definitions of compiler 
# e.g., on CC = g++-7 on OSX
ifdef PHYSICELL_CPP 
	CC := $(PHYSICELL_CPP)
endif
CC := clang++

ARCH := native # best auto-tuning
# ARCH := core2 # a reasonably safe default for most CPUs since 2007
# ARCH := corei7
# ARCH := corei7-avx # earlier i7 
# ARCH := core-avx-i # i7 ivy bridge or newer 
# ARCH := core-avx2 # i7 with Haswell or newer
# ARCH := nehalem
# ARCH := westmere
# ARCH := sandybridge # circa 2011
# ARCH := ivybridge   # circa 2012
# ARCH := haswell     # circa 2013
# ARCH := broadwell   # circa 2014
# ARCH := skylake     # circa 2015
# ARCH := bonnell
# ARCH := silvermont
# ARCH := skylake-avx512
# ARCH := nocona #64-bit pentium 4 or later 

# CFLAGS := -march=$(ARCH) -Ofast -s -fomit-frame-pointer -mfpmath=both -fopenmp -m64 -std=c++11
#CFLAGS := -march=$(ARCH) -O3 -fomit-frame-pointer -mfpmath=both -fopenmp -m64 -std=c++11
#CFLAGS := -march=$(ARCH) -O3 -fomit-frame-pointer -mfpmath=both -fopenmp -m64 -std=c++11 -I..

CFLAGS := -march=$(ARCH) -g  -fomit-frame-pointer -Xpreprocessor -fopenmp -m64 -std=c++11 -I..
PATH_TO_OMP := -L/usr/local/opt/libomp/lib

COMPILE_COMMAND := $(CC) $(CFLAGS) 

MYPATH=../..
BioFVM_OBJECTS := $(MYPATH)/BioFVM_vector.o $(MYPATH)/BioFVM_mesh.o $(MYPATH)/BioFVM_microenvironment.o $(MYPATH)/BioFVM_solvers.o $(MYPATH)/BioFVM_matlab.o \
$(MYPATH)/BioFVM_utilities.o $(MYPATH)/BioFVM_basic_agent.o $(MYPATH)/BioFVM_MultiCellDS.o $(MYPATH)/BioFVM_agent_container.o 

PhysiCell_core_OBJECTS := $(MYPATH)/PhysiCell_phenotype.o $(MYPATH)/PhysiCell_cell_container.o $(MYPATH)/PhysiCell_standard_models.o \
  $(MYPATH)/PhysiCell_cell.o $(MYPATH)/PhysiCell_custom.o $(MYPATH)/PhysiCell_utilities.o 

PhysiCell_module_OBJECTS := $(MYPATH)/PhysiCell_SVG.o $(MYPATH)/PhysiCell_pathology.o $(MYPATH)/PhysiCell_MultiCellDS.o $(MYPATH)/PhysiCell_various_outputs.o \
$(MYPATH)/PhysiCell_pugixml.o $(MYPATH)/PhysiCell_settings.o

# put your custom objects here (they should be in the custom_modules directory)

PhysiCell_custom_module_OBJECTS := custom.o

pugixml_OBJECTS := $(MYPATH)/pugixml.o

PhysiCell_OBJECTS := $(BioFVM_OBJECTS)  $(pugixml_OBJECTS) $(PhysiCell_core_OBJECTS) $(PhysiCell_module_OBJECTS)
ALL_OBJECTS := $(PhysiCell_OBJECTS) $(PhysiCell_custom_module_OBJECTS)

# compile the project 

all: main.cpp $(ALL_OBJECTS)
#	$(COMPILE_COMMAND) -o $(PROGRAM_NAME) $(ALL_OBJECTS) main.cpp 
	$(COMPILE_COMMAND) $(PATH_TO_OMP) -lomp -o $(PROGRAM_NAME) $(ALL_OBJECTS) main.cpp 

# PhysiCell core components	

PhysiCell_phenotype.o: ./core/PhysiCell_phenotype.cpp
	$(COMPILE_COMMAND) -c ./core/PhysiCell_phenotype.cpp
	
PhysiCell_digital_cell_line.o: ./core/PhysiCell_digital_cell_line.cpp
	$(COMPILE_COMMAND) -c ./core/PhysiCell_digital_cell_line.cpp

PhysiCell_cell.o: $(MYPATH)/core/PhysiCell_cell.cpp
	$(COMPILE_COMMAND) -c $(MYPATH)/core/PhysiCell_cell.cpp 

PhysiCell_cell_container.o: $(MYPATH)/core/PhysiCell_cell_container.cpp
	$(COMPILE_COMMAND) -c $(MYPATH)/core/PhysiCell_cell_container.cpp 
	
PhysiCell_standard_models.o: $(MYPATH)/core/PhysiCell_standard_models.cpp
	$(COMPILE_COMMAND) -c $(MYPATH)/core/PhysiCell_standard_models.cpp 
	
PhysiCell_utilities.o: $(MYPATH)/core/PhysiCell_utilities.cpp
	$(COMPILE_COMMAND) -c $(MYPATH)/core/PhysiCell_utilities.cpp 
	
PhysiCell_custom.o: $(MYPATH)/core/PhysiCell_custom.cpp
	$(COMPILE_COMMAND) -c $(MYPATH)/core/PhysiCell_custom.cpp 
	
# BioFVM core components (needed by PhysiCell)
	
BioFVM_vector.o: $(MYPATH)/BioFVM/BioFVM_vector.cpp
	$(COMPILE_COMMAND) -c $(MYPATH)/BioFVM/BioFVM_vector.cpp 

BioFVM_agent_container.o: $(MYPATH)/BioFVM/BioFVM_agent_container.cpp
	$(COMPILE_COMMAND) -c $(MYPATH)/BioFVM/BioFVM_agent_container.cpp 
	
BioFVM_mesh.o: $(MYPATH)/BioFVM/BioFVM_mesh.cpp
	$(COMPILE_COMMAND) -c $(MYPATH)/BioFVM/BioFVM_mesh.cpp 

BioFVM_microenvironment.o: $(MYPATH)/BioFVM/BioFVM_microenvironment.cpp
	$(COMPILE_COMMAND) -c $(MYPATH)/BioFVM/BioFVM_microenvironment.cpp 

BioFVM_solvers.o: $(MYPATH)/BioFVM/BioFVM_solvers.cpp
	$(COMPILE_COMMAND) -c $(MYPATH)/BioFVM/BioFVM_solvers.cpp 

BioFVM_utilities.o: $(MYPATH)/BioFVM/BioFVM_utilities.cpp
	$(COMPILE_COMMAND) -c $(MYPATH)/BioFVM/BioFVM_utilities.cpp 
	
BioFVM_basic_agent.o: $(MYPATH)/BioFVM/BioFVM_basic_agent.cpp
	$(COMPILE_COMMAND) -c $(MYPATH)/BioFVM/BioFVM_basic_agent.cpp 
	
BioFVM_matlab.o: $(MYPATH)/BioFVM/BioFVM_matlab.cpp
	$(COMPILE_COMMAND) -c $(MYPATH)/BioFVM/BioFVM_matlab.cpp

BioFVM_MultiCellDS.o: $(MYPATH)/BioFVM/BioFVM_MultiCellDS.cpp
	$(COMPILE_COMMAND) -c $(MYPATH)/BioFVM/BioFVM_MultiCellDS.cpp
	
pugixml.o: $(MYPATH)/BioFVM/pugixml.cpp
	$(COMPILE_COMMAND) -c $(MYPATH)/BioFVM/pugixml.cpp
	
# standard PhysiCell modules

PhysiCell_SVG.o: $(MYPATH)/modules/PhysiCell_SVG.cpp
	$(COMPILE_COMMAND) -c $(MYPATH)/modules/PhysiCell_SVG.cpp

PhysiCell_pathology.o: $(MYPATH)/modules/PhysiCell_pathology.cpp
	$(COMPILE_COMMAND) -c $(MYPATH)/modules/PhysiCell_pathology.cpp

PhysiCell_MultiCellDS.o: $(MYPATH)/modules/PhysiCell_MultiCellDS.cpp
	$(COMPILE_COMMAND) -c $(MYPATH)/modules/PhysiCell_MultiCellDS.cpp

PhysiCell_various_outputs.o: $(MYPATH)/modules/PhysiCell_various_outputs.cpp
	$(COMPILE_COMMAND) -c $(MYPATH)/modules/PhysiCell_various_outputs.cpp

PhysiCell_pugixml.o: $(MYPATH)/modules/PhysiCell_pugixml.cpp
	$(COMPILE_COMMAND) -c $(MYPATH)/modules/PhysiCell_pugixml.cpp
	
PhysiCell_settings.o: $(MYPATH)/modules/PhysiCell_settings.cpp
	$(COMPILE_COMMAND) -c $(MYPATH)/modules/PhysiCell_settings.cpp
	
# user-defined PhysiCell modules

custom.o: ./custom_modules/custom.cpp 
	$(COMPILE_COMMAND) -c ./custom_modules/custom.cpp

# cleanup

reset:
	rm -f *.cpp 
	cp ./sample_projects/Makefile-default Makefile 
	rm -f ./custom_modules/*
	touch ./custom_modules/empty.txt 
	touch ALL_CITATIONS.txt 
	rm ALL_CITATIONS.txt 
	cp ./config/PhysiCell_settings_default.xml ./config/PhysiCell_settings.xml 
	touch ./config/empty.txt
	
clean:
	rm -f *.o
	rm -f $(PROGRAM_NAME)*
	
data-cleanup:
	rm -f *.mat
	rm -f *.xml
	rm -f *.svg
	rm -f ./output/*
	touch ./output/empty.txt
	touch ./config/empty.txt
	
# archival 
	
zip:
	zip -r latest.zip Makefile* *.cpp *.h BioFVM/* config/* core/* custom_modules/* matlab/* modules/* sample_projects/* 
	cp latest.zip $$(date +%b_%d_%Y_%H%M).zip
	cp latest.zip VERSION_$(VERSION).zip 
	mv *.zip archives/
	
tar:
	tar --ignore-failed-read -czf latest.tar Makefile* *.cpp *.h BioFVM/* config/* core/* custom_modules/* matlab/* modules/* sample_projects/* 
	cp latest.tar $$(date +%b_%d_%Y_%H%M).tar
	cp latest.tar VERSION_$(VERSION).tar
	mv *.tar archives/

unzip: 
	cp ./archives/latest.zip . 
	unzip latest.zip 
	
untar: 
	cp ./archives/latest.tar .
	tar -xzf latest.tar
