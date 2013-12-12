CORDIR = correlate
PRECOMPDIR = precompute
ESTDIR = estimators
BINDIR = bin

help:
	@echo " "
	@echo "Use one of the targets below"
	@echo "	all:		build all programs"
	@echo "	pre:		build precompute code"
	@echo "	corr:		build correlation code"
	@echo "	est:		build estimation code"
	@echo "	clean:		clean all builds"

all: corr pre est

corr:
	@cd $(CORDIR); make; mv correlate ../$(BINDIR)/;

pre:
	@cd $(PRECOMPDIR); make; mv precompute ../$(BINDIR)/;

est:
	@cd $(ESTDIR); make; mv estimate ../$(BINDIR)/;

clean: clean_corr clean_pre clean_est clean_bin

clean_est:
	@cd $(ESTDIR); make clean;

clean_corr:
	@cd $(CORDIR); make clean;

clean_pre:
	@cd $(PRECOMPDIR); make clean;

clean_bin:
	@cd $(BINDIR);rm *;
