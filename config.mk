
LPADAPT_PATH ?= $(HOME)/LpAdaptation_Code/

PSSA_PREFIX ?= $(HOME)/.local/

BOOST_PYTHON_LIBS ?= -lboost_python38 -lboost_numpy38

# MEX params
MEX_EXT := $(shell mexext)

MEX_MCC := mcc
MEX_MEX := mex

TPUT := $(shell which tput 2>/dev/null)
ifneq ($(TPUT), "")
  # COLORS
  GREEN  := $(shell tput -Txterm setaf 2)
  YELLOW := $(shell tput -Txterm setaf 3)
  WHITE  := $(shell tput -Txterm setaf 7)
  RESET  := $(shell tput -Txterm sgr0)
else
  GREEN  := 
  YELLOW := 
  WHITE  := 
  RESET  := 
endif
