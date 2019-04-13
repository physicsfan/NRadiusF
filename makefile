CC=gfortran
CFLAGS=-c -Wall -g
LDFLAGS=
SOURCES=nradius.f radial.f potentials.f headings.f getdata.f physicsl_constants.f
OBJECTS=$(SOURCES:.f=.o)
EXECUTABLE=nradius

all: $(SOURCES) $(EXECUTABLE)

$(EXECUTABLE): $(OBJECTS)
	$(CC) $(LDFLAGS) $(OBJECTS) -o $@

.f.o:
	$(CC) $(CFLAGS) $< -o $@

.PHONY: clean

clean:
# For Windows
	del $(OBJECTS)
# for Linux
#	rm $(OBJECTS)
