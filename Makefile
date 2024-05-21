INSTALLDIR = test_install

toinstall = webservice.py

default :
	@echo Do "make install INSTALLDIR=<dir>"
	@echo Dev : make install INSTALLDIR=/global/cfs/cdirs/m4616/desi-gaia-server-rknop-dev
	@echo Production : make install INSTALLDIR=/global/cfs/cdirs/m4616/desi-gaia-server

install : $(patsubst %, $(INSTALLDIR)/%, $(toinstall))

$(INSTALLDIR)/% : %
	install -Dcp $< $@
