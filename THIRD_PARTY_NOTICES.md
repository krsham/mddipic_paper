# Algorithm sources and MATLAB dependencies

## OTFS transforms

`otfs_modulate` and `otfs_demodulate` implement the rectangular-pulse inverse Zak/Zak equations: `T = ifft(XDD,[],2)/M` and `XDD = M*fft(T,[],2)`. The `helperOTFSmod` and `helperOTFSdemod` files are compatibility wrappers for these project implementations.

Mathematical reference: MathWorks, [OTFS Modulation](https://www.mathworks.com/help/comm/ug/otfs-modulation.html). Earlier research code used the MathWorks example helpers; those source files are not included here.

## Installed MATLAB functionality

| Function/facility | Use | Dependency |
|---|---|---|
| `delayseq` | Finite-frame delay in physical channels, estimator atoms and the reference delay matrix | Phased Array System Toolbox; [documentation](https://www.mathworks.com/help/phased/ref/delayseq.html) |
| `qammod`, `qamdemod` | QPSK pilots and 16-QAM payload | Communications Toolbox |
| `awgn` | Optional scalar-SNR channel mode | Communications Toolbox |
| `matlab.System` and lifecycle methods | Estimator objects | MATLAB |
| FFTs, linear algebra, `pcg`, `lsqr`, `inputParser` | Transforms, fitting and equalization | MATLAB |
| `parfor` and local pool | Optional parallelism | Parallel Computing Toolbox |

`makeDelayseqReferenceKernel_v1` obtains its matrix by calling installed `delayseq` on separate complex basis vectors. Forward and adjoint operations use that matrix and its conjugate transpose. Each construction checks forward parity and the adjoint identity. Storage is quadratic in frame length.

MATLAB binaries, toolbox implementation source and example helper copies are not distributed. The additional MATLAB permission covers combination with an authorized MATLAB installation; it does not supply a MATLAB/toolbox license.

## Published methods and channel models

The baseline files are project implementations/adaptations of the cited methods. NOMP follows the Mamandipoor–Ramasamy–Madhow construction and Li 3-D adaptation, with fractional-delay refinement for the common waveform. DDAFR/MDDAFR use the Zacharia–Vani Devi search structures on the shared sampled observation model. OMP follows the cited sparse OTFS formulation. Their interfaces, gain rules and shared model are documented in LIBRARY_STRUCTURE.md.

The TDL profile tables derive from 3GPP TR 38.901. The receive-angle construction is a project extension to those scalar delay/power profiles. Full references appear in REFERENCES.bib.
