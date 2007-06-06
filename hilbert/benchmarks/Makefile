M3 = 250000 500000 750000 1000000 1500000 2000000 2500000 3000000
M6 = $(M3) 3500000 4000000 4500000 5000000 5500000 6000000

all: nn hx ic

# store vs. forget
TSFZ = $(M3:%=z8_%.dat)
TSFF = $(M3:%=f8_%.dat)
TSFR = $(M3:%=real_%.dat)

SZ: $(TSFZ:.dat=.sz.flag)
SF: $(TSFF:.dat=.sf.flag)
FZ: $(TSFZ:.dat=.fz.flag)
FF: $(TSFF:.dat=.ff.flag)
SR: $(TSFR:.dat=.sr.flag)
FR: $(TSFR:.dat=.fr.flag)
SRM: $(TSFR:.dat=.srm.flag)
FRM: $(TSFR:.dat=.frm.flag)

sf: SZ SF FZ FF SR FR SRM FRM

%.sz.flag: %.dat
	../tree_test -s 5 -l store-z.csv -t 100 -q 100 -m 5 $*.dat $*.upd || touch $@.failed
	touch $@

%.sf.flag: %.dat
	../tree_test -s 5 -l store-f.csv -t 100 -q 100 -m 5 $*.dat $*.upd || touch $@.failed
	touch $@

%.fz.flag: %.dat
	../tree_test -s 1 -l forget-z.csv -t 100 -q 100 -m 1 -f $*.dat $*.upd || touch $@.failed
	touch $@

%.ff.flag: %.dat
	../tree_test -s 1 -l forget-f.csv -t 100 -q 100 -m 1 -f $*.dat $*.upd || touch $@.failed
	touch $@

%.sr.flag: %.dat
	../tree_test -s 5 -l store-r.csv -t 100 -q 100 -m 5 $*.dat $*.upd || touch $@.failed
	touch $@

%.fr.flag: %.dat
	../tree_test -s 1 -l forget-r.csv -t 100 -q 100 -m 1 -f $*.dat $*.upd || touch $@.failed
	touch $@

%.srm.flag: %.dat
	../tree_test -s 1 -l mem-store-r.csv -t 1 -q 1 $*.dat || touch $@.failed
	touch $@

%.frm.flag: %.dat
	../tree_test -s 1 -l mem-forget-r.csv -t 1 -q 1 -f $*.dat || touch $@.failed
	touch $@

# normalization vs. non-normalization
TNZ = $(M3:%=z6_%.dat)
TNF = $(M3:%=f6_%.dat)
TNFZ = $(M3:%=f2z4_%.dat)
TNR = $(M3:%=real_%.dat)

NNZ: $(TNZ:.dat=.nnz.flag)
NNF: $(TNF:.dat=.nnf.flag)
NF: $(TNF:.dat=.nf.flag)
NNFZ: $(TNFZ:.dat=.nnfz.flag)
NFZ: $(TNFZ:.dat=.nfz.flag)
NNR: $(TNR:.dat=.nnr.flag)
NR: $(TNR:.dat=.nr.flag)

nn: NNZ NNF NF NNFZ NFZ NNR NR

%.nnz.flag: %.dat
	../tree_test -s 10 -l nonorm-z.csv -t 1 -q 1 $*.dat || touch $@.failed
	touch $@

%.nnf.flag: %.dat
	../tree_test -s 10 -l nonorm-f.csv -t 1 -q 1 $*.dat || touch $@.failed
	touch $@

%.nnfz.flag: %.dat
	../tree_test -s 10 -l nonorm-fz.csv -t 1 -q 1 $*.dat || touch $@.failed
	touch $@

%.nnr.flag: %.dat
	../tree_test -s 10 -l nonorm-r.csv -t 1 -q 1 $*.dat || touch $@.failed
	touch $@

%.nf.flag: %.dat
	../tree_test -s 10 -l norm-f.csv -t 1 -q 1 -n $*.dat || touch $@.failed
	touch $@

%.nfz.flag: %.dat
	../tree_test -s 10 -l norm-fz.csv -t 1 -q 1 -n $*.dat || touch $@.failed
	touch $@

%.nr.flag: %.dat
	../tree_test -s 10 -l norm-r.csv -t 1 -q 1 -n $*.dat || touch $@.failed
	touch $@

CTF = $(M6:%=f4_%.dat)
CTFZ = $(M6:%=f2z2_%.dat)
CTR = $(M6:%=real4_%.dat)

IF: $(CTF:.dat=.if.flag)
CF: $(CTF:.dat=.cf.flag)
IFZ: $(CTFZ:.dat=.ifz.flag)
CFZ: $(CTFZ:.dat=.cfz.flag)
IR: $(CTR:.dat=.ir.flag)
CR: $(CTR:.dat=.cr.flag)

ic: IF CF IFZ CFZ IR CR

%.ir.flag: %.dat
	../tree_test -s 10 -l iter-r.csv -t 100 -q 100 -m 10 $*.dat $*.upd || touch $@.failed
	touch $@
	
%.cr.flag: %.dat
	../tree_test -s 10 -l const-r.csv -t 100 -q 100 -m 10 -c $*.dat $*.upd || touch $@.failed
	touch $@
	
%.if.flag: %.dat
	../tree_test -s 10 -l iter-f.csv -t 100 -q 100 -m 10 $*.dat $*.upd || touch $@.failed
	touch $@
	
%.cf.flag: %.dat
	../tree_test -s 10 -l const-f.csv -t 100 -q 100 -m 10 -c $*.dat $*.upd || touch $@.failed
	touch $@
	
%.ifz.flag: %.dat
	../tree_test -s 10 -l iter-fz.csv -t 100 -q 100 -m 10 $*.dat $*.upd || touch $@.failed
	touch $@
	
%.cfz.flag: %.dat
	../tree_test -s 10 -l const-fz.csv -t 100 -q 100 -m 10 -c $*.dat $*.upd || touch $@.failed
	touch $@

HXZ = $(M6:%=z8_%.dat)
HXF = $(M6:%=f8_%.dat)
HXFZ = $(M6:%=f2z6_%.dat)
HXR = $(M6:%=real_%.dat)
	
hz: $(HXZ:.dat=.hz.flag)
xz: $(HXZ:.dat=.xz.flag)
hf: $(HXF:.dat=.hf.flag)
xf: $(HXF:.dat=.xf.flag)
hfz: $(HXFZ:.dat=.hfz.flag)
xfz: $(HXFZ:.dat=.xfz.flag)
hr: $(HXF:.dat=.hr.flag)
xr: $(HXF:.dat=.xr.flag)

hx: hz xz hf xf hfz xfz hr xr

%.hr.flag: %.dat
	../tree_test -s 10 -l hilbert-r.csv -t 100 -q 100 -m 10 $*.dat || touch $@.failed
	touch $@

%.xr.flag: %.dat
	../tree_test -s 10 -l xyz-r.csv -t 100 -q 100 -m 10 -x $*.dat || touch $@.failed
	touch $@

%.hz.flag: %.dat
	../tree_test -s 10 -l hilbert-z.csv -t 100 -q 100 -m 10 $*.dat || touch $@.failed
	touch $@

%.xz.flag: %.dat
	../tree_test -s 10 -l xyz-z.csv -t 100 -q 100 -m 10 -x $*.dat || touch $@.failed
	touch $@

%.hf.flag: %.dat
	../tree_test -s 10 -l hilbert-f.csv -t 100 -q 100 -m 10 $*.dat || touch $@.failed
	touch $@

%.xf.flag: %.dat
	../tree_test -s 10 -l xyz-f.csv -t 100 -q 100 -m 10 -x $*.dat || touch $@.failed
	touch $@

%.hfz.flag: %.dat
	../tree_test -s 10 -l hilbert-fz.csv -t 100 -q 100 -m 10 $*.dat || touch $@.failed
	touch $@

%.xfz.flag: %.dat
	../tree_test -s 10 -l xyz-fz.csv -t 100 -q 100 -m 10 -x $*.dat || touch $@.failed
	touch $@

GRF = $(M3:%=f6_%.dat)
GRF8 = $(M3:%=f8_%.dat)
GRR = $(M3:%=real_%.dat)
	
grt: $(GRR:.dat=.gr.flag)
rrt: $(GRR:.dat=.rr.flag)

gr: rrt grt
RTREE = $(HOME)/spatialindex/regressiontest/rtree/rtree

%.rr.flag: %.dat
	$(RTREE) $*.dat 100 rtree-r.csv || touch $@.failed
	touch $@

%.gr.flag: %.dat
	../tree_test -s 1 -l geocube-r.csv -t 10 -q 100 $*.dat || touch $@.failed
	touch $@
