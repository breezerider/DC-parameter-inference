include ../config.mk

.ALL: pyPSSAlib.so

.PHONY: clean pyPSSAlib_clean

pyPSSAlib.so : pyPSSAlib.cpp.o
	g++ -fPIC -g  -shared  -o $@ $< -L$(PSSA_PREFIX)/lib/ -lpssa $(BOOST_PYTHON_LIBS) $(shell python3-config --libs)

pyPSSAlib.cpp.o : pyPSSAlib.cpp
	g++ -fPIC -g $(shell python3-config --includes) -I$(PSSA_PREFIX)/include/ -c $< -o $@

pyPSSAlib_clean:
	-rm pyPSSAlib.cpp.o pyPSSAlib.so

clean: pyPSSAlib_clean
