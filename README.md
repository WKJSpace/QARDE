# QARDE

**QARDE** is an adaptive and reconfigurable non-binary LDPC (NB-LDPC) decoder engine for **continuous-variable quantum key distribution (CV-QKD)**, targeting **RFSoC FPGA systems**. It is designed in **C/C++ High-Level Synthesis (HLS)** and integrates three representative decoding algorithms within a unified architecture:

- **EMS** (Extended Min-Sum)
- **MM** (Min-Max)
- **FB-EMS** (Forward-Backward EMS)

QARDE is designed to balance **error-correction performance**, **throughput**, and **hardware cost** across different operating conditions, especially the low-SNR regime relevant to CV-QKD.

---

## Overview

In CV-QKD systems, reconciliation is a major bottleneck, and NB-LDPC decoding dominates the computational cost at low SNR. QARDE addresses this problem with a scalable FPGA architecture that combines:

- deeply pipelined message passing
- fine-grained parallelism
- carefully organized on-chip memory
- runtime algorithm configurability

The architecture supports multiple code sizes and field orders, and is implemented on **RFSoC4×2** platforms for practical CV-QKD deployment.
---

## Key Features

- **Multi-algorithm decoder framework**  
  Unified support for **EMS**, **MM**, and **FB-EMS**.

- **Runtime configurability**  
  QARDE provides interfaces to configure:
  - GF(q) field order
  - NB-LDPC code parameters
  - maximum decoding iterations
  - algorithm mode selection

- **HLS-based FPGA design**  
  Implemented using a C/C++ HLS flow to support modular design, architectural exploration, and system integration.

- **Scalable parallelism**  
  Uses a multi-level parallelism mechanism to trade off throughput and resource usage across FPGA targets.

- **CV-QKD integration**  
  Designed as part of a practical CV-QKD reconciliation system using RFSoC DAC/ADC interfaces. The system diagram is shown in **Fig. 4 on page 9** of the paper. 

- **Open-source release**  
  The paper explicitly states that QARDE is released as an open-source framework. 

---

## Decoder Flow

QARDE follows a common NB-LDPC message-passing flow:

1. **Initialization**
2. **Variable Node Update (VNU)**
3. **Permutation**
4. **Check Node Update (CNU)**
5. **Inverse Permutation**
6. **Decision / syndrome check**

This unified flow is shared across EMS, MM, and FB-EMS, while the CNU stage changes depending on the selected decoding mode.
---

## Supported Decoder Modes

### EMS
EMS truncates each q-ary message to the top-`nm` candidates and performs reduced-complexity check-node processing. It achieves near-BP performance with significantly lower complexity.

### MM
MM uses full-length q-ary messages and replaces sum-based convolution with min-max operations, resulting in a regular and hardware-friendly implementation with high throughput.

### FB-EMS
FB-EMS keeps truncated lists like EMS, but replaces the direct multi-input EMS convolution with a forward-backward chain of 2-input EMS operations. This reduces logic and memory usage while preserving strong decoding performance.

---

## Main Architecture Modules

QARDE is organized into six main modules:

- **QARDE_init**  
  Builds the Tanner graph from COO-format `H` data, constructs adjacency tables, edge metadata, permutation tables, and clears message buffers.

- **QARDE_vnu**  
  Performs variable-node update, extrinsic combination, minimum-based normalization, optional scaling/offset correction, and LDR clipping.

- **QARDE_perm**  
  Maps VN-domain messages into CN-domain messages using edge coefficient permutation tables.

- **QARDE_cnu**  
  Implements the mode-dependent check-node update for EMS, MM, or FB-EMS.

- **QARDE_invperm**  
  Converts C2V messages back from the CN domain to the VN domain.

- **QARDE_dec**  
  Computes a-posteriori metrics, produces hard decisions, and performs syndrome checking for early stopping. 
