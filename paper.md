# GIRFReco.jl: An open-source pipeline for spiral MRI Reconstruction in Julia

## Motivation

High-quality image reconstruction from frequency-domain data obtained during an MR image acquisition requires accurate characterization of the system hardware, for example via a gradient impulse response function (GIRF) [1]. If unaccounted for, system imperfections induce image artifacts such as ghosting, blurring or geometric distortion. This is important in the case of fast, SNR-efficient imaging protocols, especially for those employing spiral readout gradient waveforms [2], which place high demands on the hardware of the MRI system [3]. 

Furthermore, fast and efficient imaging necessitates acquisition in a single or few excitations, which prolongs readouts and increases their susceptibility to static off-resonance (B0) inhomogeneity. For non-Cartesian imaging, the resulting artifacts can be accounted for by incorporating off-resonance maps into the image reconstruction [4].

Therefore, state-of-the-art spiral imaging approaches [5,6] often rely on an expanded signal model incorporating system imperfections and off-resonance maps, in combination with parallel imaging acceleration using multiple receiver coils and the respective iterative non-Cartesian image reconstruction algorithms (e.g., cg-SENSE [7]).

_GIRFReco.jl_ aims to provide an open-source, single ecosystem implementation (in Julia [8]) of this approach. In particular, to enable robust and fast MRI with spiral gradient waveforms, _GIRFReco.jl_ is designed as an end-to-end signal processing pipeline from raw open-standard MR data ([ISMR]MRD [9]) to final reconstructed images (NIfTI format commonly used in processing packages for neuroimaging). Therefore, it integrates system characterization via GIRF for accurate representation of the encoding fields with calibration data (coil sensitivity and static off-resonance maps) and non-Cartesian iterative parallel imaging reconstruction [10].

## Statement of Need

Openly available software implementations for the correction of system imperfections and static off-resonance in MRI often depend upon high-quality and mature image reconstruction suites, such as MRecon (**TODO Reference/Link to Gyrotools** ), Gadgetron [11] or BART [12]. 

However, the aforementioned complexity of the image reconstruction task for spiral MRI currently necessitates to integrate multiple tools from different software suites, developed in different programming languages (e.g., C, C++, Python, Matlab), often in a containerized form, to establish a performant and comprehensive image reconstruction workflow [13].

Extending or even understanding frameworks which are built in multiple programming languages requires simultaneous expertise across languages and environments, adding significant overhead and complexity to development. This presents a significant barrier to efficient and reproducible image reconstruction and limits software accessibility and sustainability. 

The programming language Julia [8] provides a means to solve this multiple-language problem by providing a high-level interface to low-level compiled code, i.e., enabling fast prototyping with limited resources in an academic setting, while delivering performant execution times without code rewrite, all within a single environment.

In this work, we introduce _GIRFReco.jl_ (first presented at the annual meeting of ISMRM 2021 [14]), which implements an end-to-end, self-contained processing and image reconstruction pipeline for spiral MR data completely in Julia, incorporating the different model and correction components mentioned before [4,5,6,7,10]. In the spirit of software re-usabulity and sustainbility, _GIRFReco.jl_'s core image reconstruction routines are based on MRIReco.jl [15], a comprehensive open-source image reconstruction toolbox in Julia. In particular, it utilizes MRIReco.jl's iterative non-Cartesian image reconstruction [7], sensitivity map estimation using ESPiRIT [16], off-resonance correction [17], as well as open-source handling of MR raw data ((ISMR)MRD, [9]). In the course of this work, we added robust off-resonance map estimation [18] as a component (MRIFieldmaps.jl) to MRIReco.jl.

With GIRFReco.jl, we provide a self-contained, open-source spiral image reconstruction pipeline built on MRIReco.jl and include tools which correct trajectory imperfection using a Julia implementation of the GIRF correction [1,10] (adapted from a Matlab implementation of the original authors, **TODO link to GitHub Vannesjo mrigradients**), output corrected MR raw data in ISMRMRD format [9] for the spiral acquisitions, and perform necessary intermediate processing steps. GIRFReco.jl incorporates flexible parameterization of arbitrary input k-space data and supports multi-shot acquisition. It is compatible with GIRF measurements acquired using phantom scans [19,20,21], and thus can provide model-based corrections without the need of external hardware [20]. Because of the design property of Julia to create compiled low-level code, the provided image reconstruction pipeline is fast and performant.

**TODO: Maybe paragpraph on outlook (extendability), and current supported use cases (spiral diffusion working well!)**

## Functionality

### Required Inputs

_GIRFReco.jl_ requires raw MR (k-space) data (in ISMRMRD or Siemens raw data) of the following scans as input:

1. Multi-echo Gradient-echo spin-warp (Cartesian) scan
    - at least two echo times (e.g., xx and yy (water fat in-phase at 3T))
2. Spiral scan
    - single or multi-interleaf

At the moment, the slice geometries of (1) and (2) must match.

### Overview of Components

The following components are utilized within the spiral reconstruction pipeline of _GIRFReco.jl_ (Fig. 1), and separated into different (sub-)packages.

1. Core iterative image reconstruction, using Julia package _MRIReco.jl_ 
    - CG-SENSE [7] algorithm 
    - ESPIRIT [16] for sensitivity maps
2. Model-based correction components
    - Smoothed B0 map estimation, custom implementation of [18] in _MRIFieldMaps.jl_
    - Static B0 map correction, accelerated by time-segmented implementation in _MRIReco.jl [17]
    - Gradient impulse response function (GIRF) [1]
        - measured in phantom [19]
        - open-source computation [22] 
        - prediction implemented in custom package _MRIGradients.jl_ [14]


![Figure 1](paper/GIRFReco_Components.png?raw=true "GIRFReco.jl Components")

### Detailed Processing Pipeline
_GIRFReco.jl_ executes the required steps for spiral diffusion reconstruction (shown in Figure 2) in the following order:

1.	Conversion of proprietary format, vendor-specific raw image data to an open-source raw data format ((ISMR)MRD, [9])
2.	Synchronization of the data and the k-space trajectory
3.	Model-based correction of the k-space sampling points and data using the gradient impulse response function (GIRF [1], _MRIGradients.jl_) 
4.	Coil sensitivity map estimation (ESPIRIT, [16])
5.	Off-resonance (B0) map estimation (_MRIFieldmaps.jl_, [18])
6.	Non-Cartesian, iterative parallel image reconstruction with off-resonance correction (_MRIReco.jl_, [15])

![Figure 2](paper/GIRFReco_Pipeline.png?raw=true "Template Reconstruction Pipeline")

**TODO**
- Refer to example and give a bit of detail what we see and why it's good. Ideally, the example is push button
- Functionality beyond steps: 
    - Visualization: Move to Plotly.js enables remote development (HPC) with graphical feedback (w/o X11 etc., within VSCode)
    - Reproducibility: Intermediate Outputs stored with version folder (outline our folder structure?) as NIfTI files

## Getting Started
Up-to-date information about how to install _GIRFReco.jl_, run example reconstructions and apply it to your own data can be found in the README.md provided in the GitHub repository. Further technical documentation about the API is provided at **TODO: link to documenter results**, automatically generated by _Documenter.jl_.

## Conclusion

We provide _GIRFReco.jl_, an end-to-end template for spiral image reconstruction, which is open-source from raw MR data to final image, implemented in Julia as a single development ecosystem. Following best practices of software sustainability and accessibility, we re-use code of an established MR image reconstruction package (MRIReco.jl), while extending its capability to cater for the complex use case of multiple model-based corrections necessary for high-quality spiral MRI.


## References:

1. S. J. Vannesjo et al., "Gradient system characterization by impulse response measurements with a dynamic field camera," Magn. Reson. Med., vol. 69,no. 2, pp. 583-93, Feb. 2013, doi: 10.1002/mrm.24263.
2. Y. Lee et al., “On the signal-to-noise ratio benefit of spiral acquisition in diffusion MRI,” Magnetic Resonance in Medicine, vol. 85, no. 4, pp. 1924–1937, 2021, doi: https://doi.org/10.1002/mrm.28554.
3. K. T. Block and J. Frahm, "Spiral imaging: A critical appraisal," J. Magn. Reson. Imaging, vol. 21, no. 6, pp. 657-668, Jun. 2005, doi: 10.1002/jmri.20320
4. B. P. Sutton, D. C. Noll, and J. A. Fessler, “Fast, iterative image reconstruction for MRI in the presence of field inhomogeneities,” IEEE Transactions on Medical Imaging, vol. 22, no. 2, pp. 178–188, Feb. 2003, doi: 10.1109/TMI.2002.808360.
5. B. J. Wilm, C. Barmet, M. Pavan, and K. P. Pruessmann, "Higher order reconstruction for MRI in the presence of spatiotemporal field perturbations:
Higher Order Reconstruction for MRI," Magn. Reson. Med., vol. 65, no. 6, pp. 1690-1701, Jun. 2011, doi: 10.1002/mrm.22767.
6. B. J. Wilm et al., "Diffusion MRI with concurrent magnetic field monitoring," Magn. Reson. Med., vol. 74, no. 4, pp. 925-933, 2015, doi:
10.1002/mrm.25827.
7. K. P. Pruessmann, M. Weiger, P. Börnert, and P. Boesiger, “Advances in sensitivity encoding with arbitrary k-space trajectories,” Magnetic Resonance in Medicine, vol. 46, no. 4, p. 638―651, 2001, doi: 10.1002/mrm.1241.
8. J. Bezanson, A. Edelman, S. Karpinski, and V. B. Shah, “Julia: A Fresh Approach to Numerical Computing,” SIAM Rev., vol. 59, no. 1, pp. 65–98, Jan. 2017, doi: 10.1137/141000671.
9. S. J. Inati et al., “ISMRM Raw data format: A proposed standard for MRI raw datasets,” Magnetic Resonance in Medicine, vol. 77, no. 1, pp. 411–421, 2017, doi: 10.1002/mrm.26089.
10. S. J. Vannesjo et al., "Image reconstruction using a gradient impulse response model for trajectory prediction: GIRF-Based Image Reconstruction," Magn.Reson. Med., vol. 76, no. 1, pp. 45-58, Jul. 2016, doi: 10.1002/mrm.25841.
11. M. S. Hansen and T. S. Sørensen, “Gadgetron: An open source framework for medical image reconstruction,” Magn Reson Med, vol. 69, no. 6, pp. 1768–1776, Jun. 2013, doi: 10.1002/mrm.24389.
12. M. Uecker et al., “mrirecon/bart: version 0.7.00.” Zenodo, Mar. 01, 2021. doi: 10.5281/ZENODO.592960.
13. M. Veldmann, P. Ehses, K. Chow, J.-F. Nielsen, M. Zaitsev, and T. Stöcker, “Open-source MR imaging and reconstruction workflow,” Magnetic Resonance in Medicine, vol. 88, no. 6, pp. 2395–2407, 2022, doi: 10.1002/mrm.29384.
14. A. Jaffray, Z. Wu, K. Uludağ, and L. Kasper, “Open-source model-based reconstruction in Julia: A pipeline for spiral diffusion imaging,” in Proc. Intl. Soc. Mag. Reson. Med. 30, London, England, 2022, p. 2435. [Online]. Available: https://ismrm-esmrmb-ismrt2022.us3.pathable.com/meetings/virtual/poster/fMF7RmrsyKKjWbiYN
15. T. Knopp and M. Grosser, "MRIReco./I: An MRI reconstruction tramework written in Julia," Magn. Reson. Med., vol. 86, no. 3, pp. 1,633-1646, Sep. 2021,doi:10.1002/mrm.28792.
16. M. Uecker et al., "ESPIRIT-an eigenvalue approach to autocalibrating parallel MRI: Where SENSE meets GRAPPA," Magn. Reson. Med., vol. 71, no. 3,pp. 990-1001, Mar. 2014, doi: 10.1002/mrm.24751.
17. T. Knopp, H. Eggers, H. Dahnke, J. Prestin, and J. Senegas, “Iterative Off-Resonance and Signal Decay Estimation and Correction for Multi-Echo MRI,” Medical Imaging, IEEE Transactions on, vol. 28, no. 3, pp. 394–404, 2009, doi: 10.1109/TMI.2008.2006526.
18. A. K. Funai, J. A. Fessler, D. T. B. Yeo, V. T. Olafsson, and D. C. Noll, "Regularized Field Map Estimation in MRI," IEEE Trans. Med. Imaging, vol. 27, no.10, pp. 1484-1494, Oct. 2008, doi: 10.1109/TMI.2008.923956.
19. N. N. Graedel, S. A. Hurley, S. Clare, K. L. Miller, K. P. Pressmann, and S. J. Vannesjo, "Comparison of gradient impulse response functions measured with a dynamic field camera and a phantom-based technique," Barcelona/ES, 2017, p. 378.
20. M. Stich et al., "Field camera versus phantom-based measurement of the gradient system transfer function (GSTF) with dwell time compensation,"Magn. Reson. Imaging, vol. 71, pp. 125-131, Sep. 2020, do: 10.1016/.mri.2020.06.005.
21. R. K. Robison, Z. Li, D. Wang, M. B. Ooi, and J. G. Pipe, "Correction of B0 eddy current effects in spiral MRI," Magn. Reson. Med., vol. 81, no. 4, pp.2501-2513, 2019, doi:10.1002/mrm.27583.
22. Z. Wu, A. Jaffray, J. Vannesjo, K. Uludag, and L. Kasper, “MR System Stability and Quality Control using Gradient Impulse Response Functions (GIRF),” in Proc. Intl. Soc. Mag. Reson. Med. 30, London, UK, p. 0641. [Online]. Available: https://index.mirasmart.com/ISMRM2022/PDFfiles/0641.html


9. M. Lustig, S.-J. Kim, and J. M. Pauly, "A fast method for designing time-optimal gradient waveforms for arbitrary k-space trajectories," IEEE Trans. Med. Imaging, vol. 27, no. 6, pp. 866-873, Jun. 2008, doi: 10.1109/TM1.2008.922699.

