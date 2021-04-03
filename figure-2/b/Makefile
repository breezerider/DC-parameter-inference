include ../../config.mk

CUR_PATH := $(dir $(realpath $(firstword $(MAKEFILE_LIST))))
MEX_PATH := $(CUR_PATH)/../../
MEX_MODS := mexpssa mexpssa_moments mexpssa_acf mexpssa_odes mexpssa_lyapunov mexpssa_num_timepoints
MEX_DEPS := $(foreach mod,$(MEX_MODS),$(MEX_PATH)/$(mod).$(MEX_EXT))
MEX_INCL := $(foreach inc,$(MEX_DEPS), -a $(inc))
TST_NAME := ed
BIN_PATH := $(CUR_PATH)/work/
BIN_NAME := run_$(TST_NAME).sh
OUT_PATH := $(CUR_PATH)/data/
JOB_PATH := /some/path
JOB_NAME := figure-2-c
JOB_TIME := 08:00:00
JOB_CPUS := 1
JOB_MEMS := 2048
JOB_ACCT := sample-account
JOB_MAIL := example@somewhere

.PHONY: compile run analyze clean purge help

help:
	@echo -e 'usage: ${YELLOW}make${RESET} ${GREEN}[target]${RESET}\n\nTargets:\n'
	@echo -e "$$(grep -hE '^\S+:.*##' ${MAKEFILE_LIST} | sed -e 's/:.*##\s*/:/' -e 's/^\(.\+\):\(.*\)/$$#\1##:\2/' | column -c2 -t -s : | sed -e 's/$$#/${WHITE}/' -e 's/##/${RESET}/')"
	
compile: $(BIN_PATH)/$(TST_NAME) ## compile the MATLAB script into a standalone executable

run: compile ## submit the jobs ob the cluster
	for k_s in $$(seq 101 130); do \
		for k_m in $$(echo "dc cmaes"); do \
			echo "Submitting job for sample #$${k_s} using $${k_m} ..."; \
			JOB_NAME=$(JOB_NAME)-$${k_m}-$${k_s}; \
			JOB_PATH=$(JOB_PATH)/$${JOB_NAME}; \
			mkdir -p $${JOB_PATH}; \
			cp $(BIN_PATH)/{$(TST_NAME),$(BIN_NAME)} $${JOB_PATH}; \
			sed -e "s#^exe_dir=.*#exe_dir=$${JOB_PATH}#" -i $${JOB_PATH}/$(BIN_NAME); \
			sbatch --mail-user=$(JOB_MAIL) --mail-type=ALL --time=$(JOB_TIME) --ntasks=1 --cpus-per-task=$(JOB_CPUS) --job-name="$${JOB_NAME}" --account=$(JOB_ACCT) --output="$${JOB_PATH}/stdout" --error="$${JOB_PATH}/stderr" --mem-per-cpu=$(JOB_MEMS) --hint=nomultithread --workdir="$${JOB_PATH}" $${JOB_PATH}/$(BIN_NAME) "$${EBROOTMATLAB}" nSampleID $${k_s} strMethod $${k_m} strOutputPath "$(OUT_PATH)"; \
		done; \
	done

analyze: compile ## analyze the job output
	($(BIN_PATH)/$(BIN_NAME) "$${EBROOTMATLAB}" analyze $(OUT_PATH)/*.mat)

clean: ## sanitize the build environment
	-(cd $(BIN_PATH); rm $(TST_NAME) *.log *.txt $(BIN_NAME))

purge: clean ## clean all the binaries
	-rm $(MEX_DEPS)

$(BIN_PATH)/$(TST_NAME): $(BIN_PATH)/$(TST_NAME).m ../common/*.m $(LPADAPT_PATH)/LpAdaptation/*.m $(MEX_DEPS)
	(cd $(BIN_PATH); $(MEX_MCC) -m $(TST_NAME).m -a ../../common/*.m -a $(LPADAPT_PATH)/LpAdaptation/*.m $(MEX_INCL))

%.$(MEX_EXT) : %.cpp
	$(MEX_MEX) $< -I$(PSSA_PREFIX)/include -L$(PSSA_PREFIX)/lib -lut -l:libpssa.a -lgsl -lgslcblas -v -DHAVE_CONFIG_H -g -output $@
