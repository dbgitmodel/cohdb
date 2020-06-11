#
# Licensed under the EUPL, Version 1.1 or - as soon they will be approved by
# the European Commission - subsequent versions of the EUPL (the "Licence");
# You may not use this work except in compliance with the Licence.
# You may obtain a copy of the Licence at:
#
# http://ec.europa.eu/idabc/eupl
#
# Unless required by the applicable law or agreed to in writing, software
# distributed under the Licence is distributed on an "AS IS" basis,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the Licence for the specific language governing permissions and
# limitations under the Licence.
#
# Generic Makefile for COHERENS
#
# Author - Patrick Luyten
# 
# Version - @{COHERENS}Makefile   V2.9
#
# $Date: 2018-09-28 13:01:37 +0200 (Fri, 28 Sep 2018) $
#
# $Revision: 1188 $
#

# default machine dependent parameters.
FC =
FCDEFS =
FCDEBUG =
FCIFLAGS =
FLIB_FILES =
CPP =
CPPIFLAGS =
CPPF =
CPPOPTS =
CPPDEFS =

# include file with cpp options and links for compilation
DOPTSFILE = coherensflags.cmp

# options for cpp processing and libraries
include $(DOPTSFILE)

# file paths
VPATH = .:$(SEDMOD)/source:$(BIOMOD)/source:$(PARMOD)/source:$(PHYSMOD)/source
CMPDIR = $(PHYSMOD)/comps
SCMPDIR = $(SEDMOD)/comps
BCMPDIR = $(BIOMOD)/comps
PCMPDIR = $(PARMOD)/comps

# main program 
EXEFILE = coherens
MAINFILE = Coherens_Program.o

# no target
none:
	@echo No target provided

# invoke make
include $(CMPDIR)/compilers.cmp

# make command
MAKE = make

# object files
include $(CMPDIR)/objects.cmp
include $(SCMPDIR)/objects_sed.cmp
include $(BCMPDIR)/objects_bio.cmp
include $(PCMPDIR)/objects_part.cmp
MODFILES = $(MODULES) $(MODULES_SED) $(MODULES_BIO) $(MODULES_PART)
SUBFILES = $(SUBPROGS) $(SUBPROGS_SED) $(SUBPROGS_BIO) $(SUBPROGS_PART)

# compilation flags
CPPFLAGS = $(CPPIFLAGS) $(CPPOPTS) $(CPPDEFS)
FFLAGS = $(FCIFLAGS) $(FCDEFS) $(FCDEBUG)

# clean unnecessary files
clean:
	@rm -f $(MODFILES) $(SUBFILES) $(MAINFILE) $(EXEFILE) *.mod

# create executable
$(EXEFILE): $(MODFILES) $(SUBFILES) $(MAINFILE)
	$(FC) $(FFLAGS) -o $(EXEFILE) $(MODFILES) $(SUBFILES) \
	$(MAINFILE) $(FLIB_FILES)

# dependencies
include $(CMPDIR)/dependencies.cmp
include $(SCMPDIR)/dependencies_sed.cmp
include $(BCMPDIR)/dependencies_bio.cmp
include $(PCMPDIR)/dependencies_part.cmp

# suffix rules
.SUFFIXES:
.SUFFIXES:.o .f90 .F90

.F90.f90:
	@rm -f $*.f90
	$(CPPF) $(CPPFLAGS) $< ./$*.f90
.f90.o:
	$(FC) $(FFLAGS) -c $< -o $*.o
.F90.o:
	@rm -f $*.f90
	$(CPPF) $(CPPFLAGS) $< ./$*.f90
	$(FC) $(FFLAGS) -c ./$*.f90 -o $*.o
	@rm -f $*.f90
