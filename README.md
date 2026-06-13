# QARDE

### Adaptive and Reconfigurable NB-LDPC Decoder Engine for CV-QKD on RFSoC FPGAs

<p align="center">
  <img alt="platform" src="https://img.shields.io/badge/Platform-RFSoC4x2-blue">
  <img alt="design" src="https://img.shields.io/badge/Design-C%2FC%2B%2B%20HLS-orange">
  <img alt="algorithms" src="https://img.shields.io/badge/Algorithms-EMS%20%7C%20MM%20%7C%20FB--EMS-success">
  <img alt="application" src="https://img.shields.io/badge/Application-CV--QKD-purple">
  <img alt="status" src="https://img.shields.io/badge/Open%20Source-Research%20Framework-brightgreen">
</p>

<p align="center">
  <img src="Figs/QARDEoverall.png" alt="QARDE overview" width="860">
</p>

<p align="center">
  <b>QARDE</b> is an adaptive and reconfigurable <b>non-binary LDPC (NB-LDPC) decoder engine</b>
  for <b>continuous-variable quantum key distribution (CV-QKD)</b>, targeting <b>RFSoC FPGA systems</b>.
  The current implementation provides EMS, MM, and FB-EMS decoder paths in a C/C++ HLS codebase.
</p>

---

## Highlights

- **Multi-algorithm NB-LDPC framework**  
  Supports Extended Min-Sum (EMS), Min-Max (MM), and Forward-Backward EMS (FB-EMS).

- **Runtime-compatible unified accelerator**  
  `qarde_accel` keeps the original runtime-selectable decoder interface through `qarde_cnu_mode_t`.

- **Mode-specialized synthesis tops**  
  `qarde_accel_ems`, `qarde_accel_mm`, and `qarde_accel_fbems` allow each decoder mode to be synthesized independently. This reduces synthesis scope and gives Vitis HLS more freedom to optimize each datapath.

- **FPGA-oriented datapath organization**  
  The design uses pipelined message passing, partitioned q-ary message arrays, BRAM-backed buffers, and explicit CNU implementations for each algorithm.

- **Testbench coverage for public development**  
  The repository includes C-simulation tests for the algorithm kernels and dedicated top-level smoke tests for EMS, MM, and FB-EMS.

- **CV-QKD deployment focus**  
  QARDE is designed for low-SNR information reconciliation in RFSoC-based CV-QKD systems.

---

## System Context

In CV-QKD systems, information reconciliation is a major computational bottleneck, especially in the low-SNR regime where NB-LDPC decoding dominates processing cost. QARDE addresses this with a hardware-oriented HLS implementation that combines:

- q-ary LDPC message passing
- configurable iteration control
- mode-specific check-node update logic
- parallel message processing across field symbols and graph edges
- syndrome-based early stopping

<p align="center">
  <img src="Figs/CV-QKD_Bench.png" alt="QARDE deployment in RFSoC-based CV-QKD platform" width="860">
</p>
<p align="center"><i>Deployment-oriented RFSoC-based CV-QKD system context for QARDE.</i></p>

---

## Decoder Modes

| Mode | Message Representation | Check-Node Update | Typical Advantage |
|---|---|---|---|
| **EMS** | Truncated q-ary candidate lists | Reduced EMS convolution | Strong decoding quality with lower complexity than full BP |
| **MM** | Full-length q-ary messages | Min-max update | Regular datapath and high-throughput implementation |
| **FB-EMS** | Truncated q-ary candidate lists | Forward-backward chain of 2-input EMS steps | Efficient CNU realization with reduced logic and memory pressure |

### EMS

EMS keeps only the most reliable q-ary candidates in each message and performs reduced-complexity check-node processing. In QARDE, the EMS path is implemented by `qarde_decoder_ems` and the EMS CNU in `qarde_ems_cnu.hpp`.

### MM

MM operates on full-length q-ary messages and replaces sum-product style convolution with min-max operations. In QARDE, the MM path is implemented by `qarde_decoder_mm` and the MM CNU in `qarde_mm_cnu.hpp`.

### FB-EMS

FB-EMS uses truncated lists like EMS, but decomposes high-degree check-node processing into forward and backward 2-input EMS recursions. In QARDE, the FB-EMS path is implemented by `qarde_decoder_fbems` and the FB-EMS CNU in `qarde_fbems_cnu.hpp`.

---

## Decoder Flow

<p align="center">
  <img src="Figs/decoder-flow.png" alt="Unified QARDE decoder flow" width="860">
</p>

```text
Initialization
   |
   v
Variable Node Update (VNU)
   |
   v
Permutation to check-node domain
   |
   v
Check Node Update (EMS / MM / FB-EMS)
   |
   v
Inverse permutation to variable-node domain
   |
   v
A-posteriori decision and syndrome check
```

The three decoder modes share the same graph initialization, VNU, permutation, inverse permutation, decision, and syndrome-check infrastructure. The main algorithm-specific block is the CNU stage.

---

## Current HLS Design

QARDE currently exposes four accelerator entry points:

| Top Function | Source File | Purpose |
|---|---|---|
| `qarde_accel` | `src/qarde_accel.cpp` | Runtime-selectable unified top using `qarde_cnu_mode_t` |
| `qarde_accel_ems` | `src/qarde_accel_ems.cpp` | Dedicated EMS synthesis top |
| `qarde_accel_mm` | `src/qarde_accel_mm.cpp` | Dedicated MM synthesis top |
| `qarde_accel_fbems` | `src/qarde_accel_fbems.cpp` | Dedicated FB-EMS synthesis top |

The dedicated tops use the same external data style as the unified accelerator:

- input intrinsic LLR stream: `hls::stream<LLR_TYPE>`
- output hard decisions: AXI `m_axi` pointer
- control parameters: AXI-Lite
- syndrome result: AXI-Lite output

This split keeps the runtime-compatible interface available while making mode-by-mode HLS synthesis faster and easier to inspect.

---

## Main Architecture Modules

| Module | Role |
|---|---|
| `qarde_init.hpp` | Builds Tanner-graph adjacency, edge metadata, permutation tables, and clears message buffers |
| `qarde_vnu.hpp` | Performs variable-node update, extrinsic combination, correction, normalization, and clipping |
| `qarde_perm.hpp` | Permutes VN-domain messages into CN-domain order using GF coefficient tables |
| `qarde_ems_cnu.hpp` | EMS check-node update |
| `qarde_mm_cnu.hpp` | MM check-node update |
| `qarde_fbems_cnu.hpp` | FB-EMS forward-backward check-node update |
| `qarde_dec.hpp` | A-posteriori decision and syndrome checking |
| `qarde_gf.hpp` | GF arithmetic and coefficient permutation helpers |
| `qarde_tools.hpp` | Shared sort, top-k, and utility templates |
| `qarde_config.hpp` | Public accelerator function declarations |
| `nb_ldpc.hpp` | LDPC code parameters and parity-check matrix data |

Mode-specific decoder wrappers are provided in `qarde_ems.hpp`, `qarde_mm.hpp`, and `qarde_fbems.hpp`. The original runtime-selectable decoder wrapper is kept in `qarde.hpp`.

---

## Repository Structure

```text
QARDE/
├── Figs/
│   ├── CV-QKD_Bench.png
│   ├── QARDEoverall.png
│   ├── decoder-flow.png
│   └── ...
├── src/
│   ├── nb_ldpc.hpp
│   ├── qarde.hpp
│   ├── qarde_ems.hpp
│   ├── qarde_mm.hpp
│   ├── qarde_fbems.hpp
│   ├── qarde_accel.cpp
│   ├── qarde_accel_ems.cpp
│   ├── qarde_accel_mm.cpp
│   ├── qarde_accel_fbems.cpp
│   ├── qarde_config.hpp
│   ├── qarde_dec.hpp
│   ├── qarde_ems_cnu.hpp
│   ├── qarde_mm_cnu.hpp
│   ├── qarde_fbems_cnu.hpp
│   ├── qarde_gf.hpp
│   ├── qarde_init.hpp
│   ├── qarde_perm.hpp
│   ├── qarde_tools.hpp
│   └── qarde_vnu.hpp
├── testbench/
│   ├── run_csim_ems.tcl
│   ├── run_csim_mm.tcl
│   ├── run_csim_fbems.tcl
│   ├── run_csim_top_ems.tcl
│   ├── run_csim_top_mm.tcl
│   ├── run_csim_top_fbems.tcl
│   ├── tb_ems_algorithm.cpp
│   ├── tb_mm_algorithm.cpp
│   ├── tb_fbems_algorithm.cpp
│   ├── tb_top_ems.cpp
│   ├── tb_top_mm.cpp
│   └── tb_top_fbems.cpp
├── README.md
└── .gitignore
```

---

## Running C Simulation

The testbench scripts are written for Vitis HLS. On the development server, source the Vitis 2025.2 environment before running tests:

```bash
source /home/cad/xilinx/Vivado-2025.2/2025.2/Vitis/settings64.sh
```

Run algorithm-level C simulation:

```bash
vitis-run --mode hls --tcl testbench/run_csim_ems.tcl
vitis-run --mode hls --tcl testbench/run_csim_mm.tcl
vitis-run --mode hls --tcl testbench/run_csim_fbems.tcl
```

Run dedicated top-level smoke tests:

```bash
vitis-run --mode hls --tcl testbench/run_csim_top_ems.tcl
vitis-run --mode hls --tcl testbench/run_csim_top_mm.tcl
vitis-run --mode hls --tcl testbench/run_csim_top_fbems.tcl
```

---

## Synthesizing Dedicated Tops

The recommended development flow is to synthesize one dedicated top at a time. For example, a minimal FB-EMS synthesis Tcl script is:

```tcl
open_project qarde_fbems_hls
set_top qarde_accel_fbems
add_files src/qarde_accel_fbems.cpp -cflags "-std=c++17 -I./src"
open_solution "hls" -flow_target vivado
set_part {xcu55c-fsvh2892-2L-e}
create_clock -period 2.8 -name default
csynth_design
exit
```

Equivalent top/source pairs are:

| Mode | `set_top` | Source |
|---|---|---|
| EMS | `qarde_accel_ems` | `src/qarde_accel_ems.cpp` |
| MM | `qarde_accel_mm` | `src/qarde_accel_mm.cpp` |
| FB-EMS | `qarde_accel_fbems` | `src/qarde_accel_fbems.cpp` |

The unified runtime-selectable top can still be synthesized with `set_top qarde_accel` and `src/qarde_accel.cpp`, but the dedicated tops are usually better for mode-specific optimization and warning inspection.

---

## Implementation Notes

- **Design flow:** C/C++ High-Level Synthesis with Vitis HLS
- **Validated tool version:** Vitis 2025.2 development flow
- **Target device used by current scripts:** `xcu55c-fsvh2892-2L-e`
- **Target platform class:** RFSoC FPGA systems
- **Application:** CV-QKD information reconciliation
- **Supported modes:** EMS, MM, and FB-EMS
- **Runtime controls:** correction mode, iteration limit, damping, offset, and scaling parameters
- **Compile-time parameters:** GF order, graph size, row/column degree, list size, and unroll factors

---

## Reported Performance Snapshot

The QARDE paper reports RFSoC4x2 decoding throughput in the range of **2.26-12.82 Mbps** across evaluated code sizes, field orders, and decoder modes.

To highlight evaluation trends, the following figures can be included directly in the repository documentation.

<p align="center">
  <img src="Figs/GFEval.png" alt="GF evaluation results" width="85%">
  <img src="Figs/codesizeEval.png" alt="Code size evaluation results" width="85%">
  <img src="Figs/freqEval.png" alt="Frequency evaluation results" width="85%">
</p>

---

## Paper

**QARDE: A CV-QKD-based Adaptive Reconfigurable NB-LDPC Decoder Engine for RFSoC FPGA Systems**  
Kaijie Wei, Devanshu Garg, Ryutaro Nagai, Takao Tomono, and Hideharu Amano.

This repository is the public research artifact for the submitted QARDE manuscript.

---

## Citation

```bibtex
@article{wei2026qarde,
  author = {Wei, Kaijie and Garg, Devanshu and Nagai, Ryutaro and Tomono, Takao and Amano, Hideharu},
  title  = {QARDE: A CV-QKD-based Adaptive Reconfigurable NB-LDPC Decoder Engine for RFSoC FPGA Systems},
  journal = {ACM Transactions on Reconfigurable Technology and Systems},
  year   = {2026},
  note   = {Manuscript under review}
}
```

---

## Project Link

```text
https://github.com/WKJSpace/QARDE.git
```

