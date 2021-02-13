

RC     := root-config
ifeq ($(findstring $(MAKECMDTARGET),clean),)
ifeq ($(shell which $(RC) 2>&1 | sed -ne "s@.*/$(RC)@$(RC)@p"),$(RC))
MKARCH := $(wildcard $(shell $(RC) --etcdir)/Makefile.arch)
RCONFIG := $(wildcard $(shell $(RC) --incdir)/RConfigure.h)
endif
ifneq ($(MKARCH),)
include $(MKARCH)
else
ifeq ($(ROOTSYS),)
ROOTSYS = ..
endif
include $(ROOTSYS)/etc/Makefile.arch
endif
endif

#------------------------------------------------------------------------------

EVENTO        = Event.$(ObjSuf) EventDict.$(ObjSuf)
EVENTS        = Event.$(SrcSuf) EventDict.$(SrcSuf)
EVENTSO       = libEvent.$(DllSuf)
ifeq ($(PLATFORM),win32)
EVENTLIB      = libEvent.lib
else
EVENTLIB      = $(shell pwd)/$(EVENTSO)
endif


CREATETREEO    = create_tree.$(ObjSuf)
CREATETREES    = create_tree.$(SrcSuf)

PLOTTREEO    = plot_tree.$(ObjSuf)
PLOTTREES    = plot_tree.$(SrcSuf)

DYNAMICRANGEO    = dynamic_range.$(ObjSuf)
DYNAMICRANGES    = dynamic_range.$(SrcSuf)

EXECUTABLECREATETREE         = create_tree$(ExeSuf)
EXECUTABLEPLOTTREE         = plot_tree$(ExeSuf)
EXECUTABLEDYNAMICRANGE         = dynamic_range$(ExeSuf)

OBJS          = $(EVENTO) $(CREATETREEO) $(PLOTTREEO) $(DYNAMICRANGEO) 

PROGRAMS      = $(EXECUTABLECREATETREE) $(EXECUTABLEPLOTTREE) $(EXECUTABLEDYNAMICRANGE) 

ifeq ($(ARCH),aix5)
MAKESHARED    = /usr/vacpp/bin/makeC++SharedLib
endif

#------------------------------------------------------------------------------


all:            $(PROGRAMS)

$(EVENTSO):     $(EVENTO)
ifeq ($(ARCH),aix5)
		$(MAKESHARED) $(OutPutOpt) $@ $(LIBS) -p 0 $^
else
ifeq ($(PLATFORM),macosx)
# We need to make both the .dylib and the .so
		$(LD) $(SOFLAGS)$@ $(LDFLAGS) $^ $(OutPutOpt) $@ $(EXPLLINKLIBS)
ifneq ($(subst $(MACOSX_MINOR),,1234),1234)
ifeq ($(MACOSX_MINOR),4)
		ln -sf $@ $(subst .$(DllSuf),.so,$@)
endif
endif
else
ifeq ($(PLATFORM),win32)
		bindexplib $* $^ > $*.def
		lib -nologo -MACHINE:IX86 $^ -def:$*.def \
		   $(OutPutOpt)$(EVENTLIB)
		$(LD) $(SOFLAGS) $(LDFLAGS) $^ $*.exp $(LIBS) \
		   $(OutPutOpt)$@
else
		$(LD) $(SOFLAGS) $(LDFLAGS) $^ $(OutPutOpt) $@ $(EXPLLINKLIBS)
endif
endif
endif
		@echo "$@ done"

$(EXECUTABLECREATETREE):       $(EVENTSO) $(CREATETREEO)
		$(LD) $(LDFLAGS) -L${WCSIM_BASE_DIR} -I${WCSIM_BASE_DIR} -lWCSimRoot $(CREATETREEO) $(EVENTO) $(LIBS) $(OutPutOpt)$@
		@echo "$@ done"

$(EXECUTABLEPLOTTREE):       $(EVENTSO) $(PLOTTREEO)
		$(LD) $(LDFLAGS) $(PLOTTREEO) $(EVENTO) $(LIBS) $(OutPutOpt)$@
		@echo "$@ done"


$(EXECUTABLEDYNAMICRANGE):       $(DYNAMICRANGEO)
		$(LD) $(LDFLAGS) $(DYNAMICRANGEO) $(LIBS) $(OutPutOpt)$@
		@echo "$@ done"


clean:
		@rm -f $(EXECUTABLECREATETREE) $(EXECUTABLEPLOTTREE) $(EXECUTABLEDYNAMICRANGE) $(PROGRAMS)  $(OBJS) core *Dict.*

distclean:      clean
		@rm -f $(EXECUTABLECREATETREE) $(EXECUTABLEPLOTTREE) $(EXECUTABLEDYNAMICRANGE) $(PROGRAMS) $(EVENTSO) $(EVENTLIB) *Dict.* *.def *.exp \
		   *.root *.ps *.so *.lib *.dll *.d *.log .def \
		   *.pcm

.SUFFIXES: .$(SrcSuf)

###
 
Event.$(ObjSuf): Event.h
create_tree.$(ObjSuf): Event.h
plot_tree.$(ObjSuf): Event.h
dynamic_range.$(ObjSuf): Event.h

EventDict.$(SrcSuf): Event.h EventLinkDef.h
	@echo "Generating dictionary $@..."
#	$(ROOTCLING) -f $@ -c $^
	rootcint -f $@ -c $^

.$(SrcSuf).$(ObjSuf):
	$(CXX)  $(CXXFLAGS) -c $<

