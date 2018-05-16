FILE=config
PYTHON=`cat $(FILE) | grep 'python' | cut -d ' ' -f 3`
RSCRIPT=`cat $(FILE) | grep 'Rscript' | cut -d ' ' -f 3`

.PONEY: all

all:
	$(PYTHON) count.py
	$(RSCRIPT) plot.R