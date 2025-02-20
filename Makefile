CLA=clang++
CXX=g++
CXXFLAGS=-std=c++11 -Ofast -DFINAL_CHECK -DSPECIAL_HP -fpermissive
DEPS=src/beam_cky_parser.cc src/beam_cky_parser.h src/backtrace_iter.cc src/Utils/reader.h src/Utils/network.h src/Utils/codon.h src/Utils/utility_v.h src/Utils/common.h src/Utils/base.h 
BIN=bin/LinearDesign_2D
BIN_LEGACY=bin/LinearDesign_2D_legacy
UNAME_S := $(shell uname -s)
UNAME_M := $(shell uname -m)

# **Run both targets when `make` is executed without arguments**
all: lineardesign_2D lineardesign_2D_legacy

# Original target for normal LinearDesign
lineardesign_2D: $(DEPS)
	@echo "Compiling $(BIN) from src/linear_design.cpp ..."
	chmod +x lineardesign
	mkdir -p ./bin
	export LD_LIBRARY_PATH=.:$$LD_LIBRARY_PATH

ifeq ($(UNAME_S),Linux)
	@if $(CXX) $(CXXFLAGS) src/linear_design.cpp -o $(BIN) src/Utils/libraries/LinearDesign_linux64.so; then \
		echo "Linux system; compiled with g++; finished."; \
		echo "Compilation Succeed!"; \
	else \
		echo "Try another .so file."; \
		if $(CXX) $(CXXFLAGS) src/linear_design.cpp -o $(BIN) src/Utils/libraries/LinearDesign_linux64_old.so; then \
			echo "Linux system; compiled with g++; finished."; \
			echo "Compilation Succeed!"; \
		else \
			echo "Compilation failed! Make sure it is either Linux-64 or Mac."; \
		fi \
	fi 
else
	@if [ "$(UNAME_M)" = 'arm64' ]; then \
		if 	$(CLA) $(CXXFLAGS) src/linear_design.cpp -o $(BIN) src/Utils/libraries/LinearDesign_Mac_M1.so; then \
			echo "Mac M1 system; compiled with clang++; finished."; \
			echo "Compilation Succeed!"; \
			echo "You may encounter a pop-up message at the first run. If so, please go to System Preferences -> Security & Privacy -> General to allow LinearDesign_Mac_M1.so to open. See README.md for details."; \
		else \
			echo "Compilation failed! Make sure it is either Linux-64 or Mac."; \
		fi \
	else \
		if 	$(CLA) $(CXXFLAGS) src/linear_design.cpp -o $(BIN) src/Utils/libraries/LinearDesign_Mac_x86.so; then \
			echo "Mac x86_64 system; compiled with clang++; finished."; \
			echo "Compilation Succeed!"; \
			echo "You may encounter a pop-up message at the first run. If so, please go to System Preferences -> Security & Privacy -> General to allow LinearDesign_Mac_x86.so to open. See README.md for details."; \
		else \
			echo "Compilation failed! Make sure it is either Linux-64 or Mac."; \
		fi \
	fi
endif

# New target for LinearDesign Legacy version
lineardesign_2D_legacy: $(DEPS)
	@echo "Compiling $(BIN_LEGACY) from src/legacies/linear_design.cpp ..."
	chmod +x lineardesign
	mkdir -p ./bin
	export LD_LIBRARY_PATH=.:$$LD_LIBRARY_PATH

ifeq ($(UNAME_S),Linux)
	@if $(CXX) $(CXXFLAGS) src/legacies/linear_design.cpp -o $(BIN_LEGACY) src/Utils/libraries/LinearDesign_linux64.so; then \
		echo "Linux system; compiled with g++; finished."; \
		echo "Compilation Succeed!"; \
	else \
		echo "Try another .so file."; \
		if $(CXX) $(CXXFLAGS) src/legacies/linear_design.cpp -o $(BIN_LEGACY) src/Utils/libraries/LinearDesign_linux64_old.so; then \
			echo "Linux system; compiled with g++; finished."; \
			echo "Compilation Succeed!"; \
		else \
			echo "Compilation failed! Make sure it is either Linux-64 or Mac."; \
		fi \
	fi 
else
	@if [ "$(UNAME_M)" = 'arm64' ]; then \
		if 	$(CLA) $(CXXFLAGS) src/legacies/linear_design.cpp -o $(BIN_LEGACY) src/Utils/libraries/LinearDesign_Mac_M1.so; then \
			echo "Mac M1 system; compiled with clang++; finished."; \
			echo "Compilation Succeed!"; \
			echo "You may encounter a pop-up message at the first run. If so, please go to System Preferences -> Security & Privacy -> General to allow LinearDesign_Mac_M1.so to open. See README.md for details."; \
		else \
			echo "Compilation failed! Make sure it is either Linux-64 or Mac."; \
		fi \
	else \
		if 	$(CLA) $(CXXFLAGS) src/legacies/linear_design.cpp -o $(BIN_LEGACY) src/Utils/libraries/LinearDesign_Mac_x86.so; then \
			echo "Mac x86_64 system; compiled with clang++; finished."; \
			echo "Compilation Succeed!"; \
			echo "You may encounter a pop-up message at the first run. If so, please go to System Preferences -> Security & Privacy -> General to allow LinearDesign_Mac_x86.so to open. See README.md for details."; \
		else \
			echo "Compilation failed! Make sure it is either Linux-64 or Mac."; \
		fi \
	fi
endif

.PHONY: clean

clean:
	rm -f $(BIN) $(BIN_LEGACY)
