include ../config.mk

.ALL: fig1

.PHONY: fig1 clean

fig1: run_fig1.sh

run_fig1.sh: run_fig1.sh.in
	sed -e "s#@PREFIX@#$(PSSA_PREFIX)#g" -e "s#@SHELL@#$(SHELL)#g" $< >$@

clean:
	-rm run_fig1.sh
