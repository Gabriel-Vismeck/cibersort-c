# =============================================================================
#  Makefile do CIBERSORT em C
#
#  make              compila (build/cibersort)
#  make OPENMP=1     compila usando todos os núcleos (OpenMP)
#  make limpar       apaga a pasta build/
# =============================================================================

CC      ?= gcc
CFLAGS  ?= -O2
CFLAGS  += -std=c99 -Wall -Wextra -Iinclude
LDLIBS  += -lm

ifeq ($(OPENMP),1)
CFLAGS += -fopenmp
endif

FONTES  := $(wildcard src/*.c)
OBJETOS := $(patsubst src/%.c,build/%.o,$(FONTES))
ALVO    := build/cibersort

.PHONY: todos teste limpar

todos: $(ALVO)

$(ALVO): $(OBJETOS)
	$(CC) $(CFLAGS) -o $@ $^ $(LDLIBS)

build/%.o: src/%.c | build
	$(CC) $(CFLAGS) -MMD -MP -c -o $@ $<

build:
	mkdir -p build

limpar:
	rm -rf build

-include $(OBJETOS:.o=.d)
