CCFLAGS	+=	-fPIC

CSRCS		=	$(wildcard *.c)
OBJS		=	$(CSRCS:.c=.o)
OBJSS		=	$(CSRCS:.c=.os)
CLHD		=	$(CSRCS:.c=.h)
ODEP		=	$(CSRCS:.c=.d)

all: shared static
shared: librngstream.so
static: librngstream.a

librngstream.a: $(OBJSS)
	rm -rf $@ && ar rcs $@ $^

librngstream.so: $(OBJS)
	$(CC) $(CFLAGS) -shared -o $@ $^

clean:
	rm -f *.d
	rm -f *.o
	rm -f *.os
	rm -f *.ipos

clear: clean
	rm -f librngstream.so librngstream.a

$(ODEP): %.d: %.c %.h
	@echo "Generating dependency file $@"
	@set -e; rm -f $@
	@$(CC) -M $(CFLAGS) $< > $@.tmp
	@sed 's,\($*\)\.o[ :]*,\1.o $@ : ,g' < $@.tmp > $@
	@rm -f $@.tmp

include $(ODEP)

$(OBJS): %.o: %.c
	 $(CC) -c -o $@ -fPIC $(CFLAGS) $<

$(OBJSS): %.os: %.c
	 $(CC) -c -o $@ $(CFLAGS) $<

arv.ipos: $(OBJSS)
	$(CC) -ipo_c -o $@ $(CFLAGS) $^
