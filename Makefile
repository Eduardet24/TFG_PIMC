# Unified GNU Makefile for Grigori_codigo (QMC / PIMC / Worm)
#
# Configurable variables (without editing main.h):
#   DIM          = 1D | 2D | 3D        (default: 2D)
#   WORM         = 1 | 0               (default: 1)
#   NORMAL_MODES = 1 | 0               (default: 1)
#   TRAP         = 1 | 0               (default: 1)
#   OPT          = -O3                 (default: -O3)
#   DEBUG        = 0 | 1               (default: 0)
#   ASAN         = 0 | 1               (default: 0)
#   OPENMP       = 1 | 0               (default: 1)
#   TARGET       = MC                  (output binary path)
#   EXTRA_CPPFLAGS =                   (additional -D... flags)
#
# Examples:
#   make                               # builds ./MC (2D + Worm + NormalModes)
#   make DIM=1D WORM=0                 # builds ./MC (1D without Worm)
#   make MC_1D                         # builds ./MC_1D (1D without Worm)
#   make MC_2D_worm                    # builds ./MC_2D_worm (2D with Worm)
#   make MC_noworm                     # builds ./MC_noworm (2D without Worm)
#   make MC_3D                         # builds ./MC_3D (3D with Worm)
#   make MC_asan                       # builds ./MC_asan with AddressSanitizer

CC ?= gcc
DIM ?= 2D
WORM ?= 1
NORMAL_MODES ?= 1
TRAP ?= 1
OPT ?= -O3
DEBUG ?= 0
ASAN ?= 0
OPENMP ?= 1
TARGET ?= MC
EXTRA_CPPFLAGS ?=

EMPTY :=
SPACE := $(EMPTY) $(EMPTY)

ifeq ($(DIM),1D)
  DIM_FLAG := -DTRIAL_1D
else ifeq ($(DIM),2D)
  DIM_FLAG := -DTRIAL_2D
else ifeq ($(DIM),3D)
  DIM_FLAG := -DTRIAL_3D
else
  $(error Invalid DIM='$(DIM)'; expected 1D, 2D, or 3D)
endif

ifeq ($(WORM),1)
  WORM_FLAG := -DWORM_ALGORITHM
else ifeq ($(WORM),0)
  WORM_FLAG := -DNO_WORM_ALGORITHM
else
  $(error Invalid WORM='$(WORM)'; expected 1 or 0)
endif

ifeq ($(NORMAL_MODES),1)
  NM_FLAG := -DPIMC_HARMONIC_NORMAL_MODES
else ifeq ($(NORMAL_MODES),0)
  NM_FLAG := -DNO_PIMC_HARMONIC_NORMAL_MODES
else
  $(error Invalid NORMAL_MODES='$(NORMAL_MODES)'; expected 1 or 0)
endif

ifeq ($(TRAP),1)
  TRAP_FLAG := -DTRAP_POTENTIAL
else ifeq ($(TRAP),0)
  TRAP_FLAG := -DNO_TRAP_POTENTIAL
else
  $(error Invalid TRAP='$(TRAP)'; expected 1 or 0)
endif

OMP_FLAG := $(if $(filter 1,$(OPENMP)),-fopenmp,)

ifeq ($(ASAN),1)
  BASE_OPT := -O1 -g -fsanitize=address -fno-omit-frame-pointer
  ASAN_LD := -fsanitize=address
else ifeq ($(DEBUG),1)
  BASE_OPT := -O0 -g
  ASAN_LD :=
else
  BASE_OPT := $(OPT)
  ASAN_LD :=
endif

CPPFLAGS_EFF := $(DIM_FLAG) $(WORM_FLAG) $(NM_FLAG) $(TRAP_FLAG) $(EXTRA_CPPFLAGS)
CFLAGS_EFF := $(BASE_OPT) $(OMP_FLAG) -std=c99 -Wall -Wextra
LDFLAGS_EFF := $(OMP_FLAG) $(ASAN_LD)
LDLIBS_EFF := -lm

OPT_TAG := $(subst $(SPACE),_,$(strip $(BASE_OPT)))
CFG_TAG := $(DIM)_w$(WORM)_nm$(NORMAL_MODES)_t$(TRAP)_omp$(OPENMP)_dbg$(DEBUG)_asan$(ASAN)_opt$(OPT_TAG)_$(notdir $(CC))
ifneq ($(strip $(EXTRA_CPPFLAGS)),)
  CFG_HASH := $(shell printf '%s' '$(strip $(EXTRA_CPPFLAGS))' | cksum | awk '{print $$1}')
  CFG_TAG := $(CFG_TAG)_x$(CFG_HASH)
endif

BUILD_ROOT := build
BUILD_DIR := $(BUILD_ROOT)/$(CFG_TAG)
STAMP_DIR := $(BUILD_ROOT)/stamps
INTERNAL_BIN := $(BUILD_DIR)/MC.bin
TARGET_HASH := $(shell printf '%s' '$(abspath $(TARGET))' | cksum | awk '{print $$1}')
TARGET_STAMP := $(STAMP_DIR)/$(notdir $(TARGET))_$(TARGET_HASH).cfg

SRCS := $(sort $(filter-out tnk-igor.c,$(wildcard *.c)))
OBJS := $(patsubst %.c,$(BUILD_DIR)/%.o,$(SRCS))
DEPS := $(OBJS:.o=.d)

# Force relink/copy if the requested TARGET was last built with a different CFG_TAG
PREV_TARGET_CFG := $(if $(wildcard $(TARGET_STAMP)),$(strip $(file < $(TARGET_STAMP))),)
ifneq ($(PREV_TARGET_CFG),$(CFG_TAG))
  .PHONY: $(TARGET)
endif

.PHONY: all clean distclean help MC_1D MC_2D_worm MC_noworm MC_3D MC_asan

all: $(TARGET)

$(TARGET): $(INTERNAL_BIN)
	@mkdir -p $(dir $@) $(STAMP_DIR)
	@cp -f $< $@
	@printf '%s\n' '$(CFG_TAG)' > $(TARGET_STAMP)
	@echo "[MAKE] Built $@ ($(CFG_TAG))"

$(INTERNAL_BIN): $(OBJS)
	@mkdir -p $(dir $@)
	$(CC) $(OBJS) $(LDFLAGS_EFF) $(LDLIBS_EFF) -o $@

$(BUILD_DIR)/%.o: %.c
	@mkdir -p $(dir $@)
	$(CC) $(CPPFLAGS_EFF) $(CFLAGS_EFF) -MMD -MP -MF $(@:.o=.d) -c $< -o $@

MC_1D:
	@$(MAKE) --no-print-directory TARGET=MC_1D DIM=1D WORM=0 NORMAL_MODES=1

MC_2D_worm:
	@$(MAKE) --no-print-directory TARGET=MC_2D_worm DIM=2D WORM=1 NORMAL_MODES=1

MC_noworm:
	@$(MAKE) --no-print-directory TARGET=MC_noworm DIM=2D WORM=0 NORMAL_MODES=1

MC_3D:
	@$(MAKE) --no-print-directory TARGET=MC_3D DIM=3D WORM=1 NORMAL_MODES=1

MC_asan:
	@$(MAKE) --no-print-directory TARGET=MC_asan ASAN=1

clean:
	rm -rf $(BUILD_ROOT)

distclean: clean
	rm -f MC MC_1D MC_2D_worm MC_noworm MC_3D MC_asan

help:
	@echo "Targets: all (MC), MC_1D, MC_2D_worm, MC_noworm, MC_3D, MC_asan, clean, distclean"
	@echo "Options: DIM=1D|2D|3D WORM=1|0 NORMAL_MODES=1|0 TRAP=1|0 OPT=-O3 DEBUG=0|1 ASAN=0|1 OPENMP=1|0 TARGET=MC"

-include $(DEPS)
