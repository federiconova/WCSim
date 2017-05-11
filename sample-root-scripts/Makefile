

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


READTREEO    = read_tree.$(ObjSuf)
READTREES    = read_tree.$(SrcSuf)

CONVERTTOMASKO    = convert_to_mask.$(ObjSuf)
CONVERTTOMASKS    = convert_to_mask.$(SrcSuf)

SELECTTREEO    = select_tree.$(ObjSuf)
SELECTTREES    = select_tree.$(SrcSuf)

CORRECTENERGYO    = correct_energy.$(ObjSuf)
CORRECTENERGYS    = correct_energy.$(SrcSuf)

CORRECTENERGYCOSO    = correct_energy_cos.$(ObjSuf)
CORRECTENERGYCOSS    = correct_energy_cos.$(SrcSuf)

SELECTAFTERENERGYCORRECTIONO    = select_after_energy_correction.$(ObjSuf)
SELECTAFTERENERGYCORRECTIONS    = select_after_energy_correction.$(SrcSuf)

SELECTAFTERENERGYCOSCORRECTIONO    = select_after_energy_cos_correction.$(ObjSuf)
SELECTAFTERENERGYCOSCORRECTIONS    = select_after_energy_cos_correction.$(SrcSuf)

PLOTCONTAINEDMUONSO    = plot_contained_muons.$(ObjSuf)
PLOTCONTAINEDMUONSS    = plot_contained_muons.$(SrcSuf)

MUONHITTIMESO    = muon_hit_times.$(ObjSuf)
MUONHITTIMESS    = muon_hit_times.$(SrcSuf)

CREATETREEFROMTXTO    = create_tree_from_txt.$(ObjSuf)
CREATETREEFROMTXTS    = create_tree_from_txt.$(SrcSuf)

PLOTTREEO    = plot_tree.$(ObjSuf)
PLOTTREES    = plot_tree.$(SrcSuf)

WRITEALLHITSO    = write_all_hits.$(ObjSuf)
WRITEALLHITSS    = write_all_hits.$(SrcSuf)

WRITEALLPMTSO    = write_all_pmts.$(ObjSuf)
WRITEALLPMTSS    = write_all_pmts.$(SrcSuf)

WRITEDETECTORO    = write_detector.$(ObjSuf)
WRITEDETECTORS    = write_detector.$(SrcSuf)

EXECUTABLEREADTREE         = read_tree$(ExeSuf)
EXECUTABLECONVERTTOMASK         = convert_to_mask$(ExeSuf)
EXECUTABLESELECTTREE         = select_tree$(ExeSuf)
EXECUTABLECORRECTENERGY         = correct_energy$(ExeSuf)
EXECUTABLECORRECTENERGYCOS         = correct_energy_cos$(ExeSuf)
EXECUTABLESELECTAFTERENERGYCORRECTION         = select_after_energy_correction$(ExeSuf)
EXECUTABLESELECTAFTERENERGYCOSCORRECTION         = select_after_energy_cos_correction$(ExeSuf)
EXECUTABLEPLOTCONTAINEDMUONS         = plot_contained_muons$(ExeSuf)
EXECUTABLEMUONHITTIMES         = muon_hit_times$(ExeSuf)
EXECUTABLECREATETREEFROMTXT         = create_tree_from_txt$(ExeSuf)
EXECUTABLEPLOTTREE         = plot_tree$(ExeSuf)
EXECUTABLEWRITEALLHITS         = write_all_hits$(ExeSuf)
EXECUTABLEWRITEALLPMTS         = write_all_pmts$(ExeSuf)
EXECUTABLEWRITEDETECTOR         = write_detector$(ExeSuf)

OBJS          = $(EVENTO) $(READTREEO) $(CONVERTTOMASKO) $(SELECTTREEO) $(CORRECTENERGYO) $(CORRECTENERGYCOSO) $(SELECTAFTERENERGYCORRECTIONO) $(SELECTAFTERENERGYCOSCORRECTIONO) $(PLOTCONTAINEDMUONSO) $(MUONHITTIMESO) $(CREATETREEFROMTXTO) $(PLOTTREEO) $(WRITEALLHITSO) $(WRITEALLPMTSO) $(WRITEDETECTORO)

PROGRAMS      = $(EXECUTABLEREADTREE) $(EXECUTABLECONVERTTOMASK) $(EXECUTABLESELECTTREE) $(EXECUTABLECORRECTENERGY) $(EXECUTABLECORRECTENERGYCOS) $(EXECUTABLESELECTAFTERENERGYCORRECTION) $(EXECUTABLESELECTAFTERENERGYCOSCORRECTION) $(EXECUTABLEPLOTCONTAINEDMUONS) $(EXECUTABLEMUONHITTIMES) $(EXECUTABLECREATETREEFROMTXT) $(EXECUTABLEPLOTTREE) $(EXECUTABLEWRITEALLHITS)  $(EXECUTABLEWRITEALLPMTS)  $(EXECUTABLEWRITEDETECTOR) 

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

$(EXECUTABLEREADTREE):       $(EVENTSO) $(READTREEO)
		$(LD) $(LDFLAGS) $(READTREEO) $(EVENTO) $(LIBS) $(OutPutOpt)$@
		@echo "$@ done"

$(EXECUTABLECONVERTTOMASK):       $(EVENTSO) $(CONVERTTOMASKO)
		$(LD) $(LDFLAGS) $(CONVERTTOMASKO) $(EVENTO) $(LIBS) $(OutPutOpt)$@
		@echo "$@ done"

$(EXECUTABLESELECTTREE):       $(EVENTSO) $(SELECTTREEO)
		$(LD) $(LDFLAGS) $(SELECTTREEO) $(EVENTO) $(LIBS) $(OutPutOpt)$@
		@echo "$@ done"

$(EXECUTABLECORRECTENERGY):       $(EVENTSO) $(CORRECTENERGYO)
		$(LD) $(LDFLAGS) $(CORRECTENERGYO) $(EVENTO) $(LIBS) $(OutPutOpt)$@
		@echo "$@ done"

$(EXECUTABLECORRECTENERGYCOS):       $(EVENTSO) $(CORRECTENERGYCOSO)
		$(LD) $(LDFLAGS) $(CORRECTENERGYCOSO) $(EVENTO) $(LIBS) $(OutPutOpt)$@
		@echo "$@ done"

$(EXECUTABLESELECTAFTERENERGYCORRECTION):       $(EVENTSO) $(SELECTAFTERENERGYCORRECTIONO)
		$(LD) $(LDFLAGS) $(SELECTAFTERENERGYCORRECTIONO) $(EVENTO) $(LIBS) $(OutPutOpt)$@
		@echo "$@ done"

$(EXECUTABLESELECTAFTERENERGYCOSCORRECTION):       $(EVENTSO) $(SELECTAFTERENERGYCOSCORRECTIONO)
		$(LD) $(LDFLAGS) $(SELECTAFTERENERGYCOSCORRECTIONO) $(EVENTO) $(LIBS) $(OutPutOpt)$@
		@echo "$@ done"

$(EXECUTABLEPLOTCONTAINEDMUONS):       $(EVENTSO) $(PLOTCONTAINEDMUONSO)
		$(LD) $(LDFLAGS) $(PLOTCONTAINEDMUONSO) $(EVENTO) $(LIBS) $(OutPutOpt)$@
		@echo "$@ done"

$(EXECUTABLEMUONHITTIMES):       $(EVENTSO) $(MUONHITTIMESO)
		$(LD) $(LDFLAGS) $(MUONHITTIMESO) $(EVENTO) $(LIBS) $(OutPutOpt)$@
		@echo "$@ done"

$(EXECUTABLECREATETREEFROMTXT):       $(EVENTSO) $(CREATETREEFROMTXTO)
		$(LD) $(LDFLAGS) $(CREATETREEFROMTXTO) $(EVENTO) $(LIBS) $(OutPutOpt)$@
		@echo "$@ done"

$(EXECUTABLEPLOTTREE):       $(EVENTSO) $(PLOTTREEO)
		$(LD) $(LDFLAGS) $(PLOTTREEO) $(EVENTO) $(LIBS) $(OutPutOpt)$@
		@echo "$@ done"

$(EXECUTABLEWRITEALLHITS):       $(EVENTSO) $(WRITEALLHITSO)
		$(LD) $(LDFLAGS) $(WRITEALLHITSO) $(EVENTO) $(LIBS) $(OutPutOpt)$@
		@echo "$@ done"

$(EXECUTABLEWRITEALLPMTS):       $(EVENTSO) $(WRITEALLPMTSO)
		$(LD) $(LDFLAGS) $(WRITEALLPMTSO) $(EVENTO) $(LIBS) $(OutPutOpt)$@
		@echo "$@ done"


$(EXECUTABLEWRITEDETECTOR):       $(EVENTSO) $(WRITEDETECTORO)
		$(LD) $(LDFLAGS) $(WRITEDETECTORO) $(EVENTO) $(LIBS) $(OutPutOpt)$@
		@echo "$@ done"



clean:
		@rm -f $(EXECUTABLEREADTREE) $(EXECUTABLECONVERTTOMASK) $(EXECUTABLESELECTTREE) $(EXECUTABLECORRECTENERGY) $(EXECUTABLECORRECTENERGYCOS) $(EXECUTABLESELECTAFTERENERGYCORRECTION) $(EXECUTABLESELECTAFTERENERGYCOSCORRECTION) $(EXECUTABLEPLOTCONTAINEDMUONS) $(EXECUTABLEMUONHITTIMES) $(EXECUTABLECREATETREEFROMTXT) $(EXECUTABLEPLOTTREE) $(EXECUTABLEWRITEALLHITS) $(EXECUTABLEWRITEALLPMTS) $(EXECUTABLEWRITEDETECTOR) $(PROGRAMS)  $(OBJS) core *Dict.*

distclean:      clean
		@rm -f $(EXECUTABLEREADTREE) $(EXECUTABLECONVERTTOMASK) $(EXECUTABLESELECTTREE) $(EXECUTABLECORRECTENERGY) $(EXECUTABLECORRECTENERGYCOS) $(EXECUTABLESELECTAFTERENERGYCORRECTION) $(EXECUTABLESELECTAFTERENERGYCOSCORRECTION) $(EXECUTABLEPLOTCONTAINEDMUONS) $(EXECUTABLEMUONHITTIMES) $(EXECUTABLECREATETREEFROMTXT) $(EXECUTABLEPLOTTREE) $(EXECUTABLEWRITEALLHITS) $(EXECUTABLEWRITEALLPMTS) $(EXECUTABLEWRITEDETECTOR) $(PROGRAMS) $(EVENTSO) $(EVENTLIB) *Dict.* *.def *.exp \
		   *.root *.ps *.so *.lib *.dll *.d *.log .def \
		   *.pcm

.SUFFIXES: .$(SrcSuf)

###
 
Event.$(ObjSuf): Event.h
read_tree.$(ObjSuf): Event.h
convert_to_mask.$(ObjSuf): Event.h
select_tree.$(ObjSuf): Event.h
correct_energy.$(ObjSuf): Event.h
correct_energy_cos.$(ObjSuf): Event.h
select_after_energy_correction.$(ObjSuf): Event.h
select_after_energy_cos_correction.$(ObjSuf): Event.h
plot_contained_muons.$(ObjSuf): Event.h
muon_hit_times.$(ObjSuf): Event.h
create_tree_from_txt.$(ObjSuf): Event.h
plot_tree.$(ObjSuf): Event.h
write_all_hits.$(ObjSuf): Event.h
write_all_pmts.$(ObjSuf): Event.h
write_detector.$(ObjSuf): Event.h

EventDict.$(SrcSuf): Event.h EventLinkDef.h
	@echo "Generating dictionary $@..."
#	$(ROOTCLING) -f $@ -c $^
	rootcint -f $@ -c $^

.$(SrcSuf).$(ObjSuf):
	$(CXX)  $(CXXFLAGS) -c $<

