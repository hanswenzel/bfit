ALL	= fitter

fitter:
	@(cd src; make)

clean:
	@(cd src; make clean)
	@(cd config_file; make clean)
	@rm bexclfit mkbdecay

