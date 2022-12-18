MAKE = gmake

srcdir := $(subst config/../,,.)
top_srcdir := $(subst config/../,,..)
VPATH := $(subst config/../,,.)

BELLE_CONFIG_DIR = /belle/belle/b20090127_0910/src/config

PACKAGE = fullrecon
VERSION = 1.0

packagedir = .

PACKAGE_HEADER_DIR = $(PACKAGE)
export PACKAGE_HEADER_DIR
BELLE_ROUTE = src/anal/fullrecon/
export BELLE_ROUTE

SUBDIRS = doc fullrecon src
MAINDIRS = basf_if 

include $(BELLE_CONFIG_DIR)/Makefile.subdirs

