HTSPROGS=vcf2eig
BASEPROGS=eig2multifasta eig2phylip eigreduce distance
PROGS=$(BASEPROGS) $(HTSPROGS) abba-baba
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

abba-baba: $$@.o kmath.o
	$(CC) $(CFLAGS) $^ -o $@ $(LDLIBS) -lm

clean:
	rm -f $(PROGS) $(addsuffix .o,$(PROGS)) kmath.o
