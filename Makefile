HTSPROGS=vcf2eig
BASEPROGS=eig2phylip eigreduce
PROGS=$(BASEPROGS) $(HTSPROGS)
ZLIB=/opt/shared/zlib/1.2.8-gnu_4.8.0/lib
HTSDIR=../../../htslib
HTSFLAGS=-I$(HTSDIR)
HTSLIBS=-L$(HTSDIR) -L$(ZLIB) -lm -Wl,-static -lhts -lz -Wl,-Bdynamic -pthread
CFLAGS=-Wall -O2 -g $(HTSFLAGS)
LDLIBS=

all: $(PROGS)

.SECONDEXPANSION:
$(HTSPROGS): $$@.o
	$(CC) $(HTSFLAGS) $(CFLAGS) $^ -o $@ $(HTSLIBS) $(LDLIBS)

$(BASEPROGS): $$@.o
	$(CC) $(CFLAGS) $^ -o $@ $(LDLIBS)

clean:
	rm -f $(PROGS) $(addsuffix .o,$(PROGS))
