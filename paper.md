---
title: 'GIRFReco.jl: An Open-Source Pipeline for Spiral Magnetic Resonance Image (MRI) Reconstruction in Julia'
tags:
  - Julia
  - Magnetic Resonance Imaging
  - Non-Cartesian Image Reconstruction
  - Gradient Impulse Response Function (GIRF)
  - Off-resonance Correction
authors:
  - name: Alexander Jaffray
    orcid: 0000-0000-0000-0000
    equal-contrib: True
    corresponding: True
    affiliation: 1
  - name: Author2
    orcid: 0000-0000-0000-0000
    equal-contrib: True # (This is how you can denote equal contributions between multiple authors)
    corresponding: True # (This is how to denote the corresponding author)
    affiliation: 2
  - name: Author3
    equal-contrib: False # (This is how you can denote equal contributions between multiple authors)
    corresponding: False # (This is how to denote the corresponding author)
    affiliation: 2
  - name: Senior Author
    equal-contrib: False # (This is how you can denote equal contributions between multiple authors)
    corresponding: False # (This is how to denote the corresponding author)
    affiliation: "2, 3" # (Multiple affiliations must be quoted)
affiliations:
 - name: MRI Research Centre, University of British Columbia, Vancouver, Canada
   index: 1
 - name: Techna Institute, University Health Network, Ontario, Canada
   index: 2
 - name: Department of Medical Biophysics, University of Toronto, Canada
   index: 3
date: 31 January 2023
bibliography: paper/paper.bib

# Optional fields if submitting to a AAS journal too, see this blog post:
# https://blog.joss.theoj.org/2018/12/a-new-collaboration-with-aas-publishing
# aas-doi: 10.3847/xxxxx <- update this with the DOI from AAS once you know it.
# aas-journal: Astrophysical Journal <- The name of the AAS journal.
---

# Summary

Magnetic Resonance Imaging (MRI) acquires data in the frequency domain (k-space) with various possibilities of trajectories to traverse it. Traditional Cartesian MRI is notoriously slow since it's using a line-by-line pattern to cover k-space due to its intrinsically low signal-to-noise ratio (SNR). On the other hand, non-Cartesian (i.e., non-rectilinear) trajectories, such as spiral readouts, offer considerable acceleration due their high rate of measured samples per time (acquisition efficiency) [@lee_signalnoise_2021].

However, the actual k-space trajectory during MRI scanning deviates from the theoretical one due to hardware imperfections, which induces image artifacts such as ghosting, blurring or geometric distortion. This problem is exacerbated in many non-Cartesian trajectories, including spiral readouts, because these fast imaging protocols place high demands on the hardware of the MRI system [@block_spiral_2005].

In this case, an accurate characterization of the system hardware for k-space trajectory correction, for example, via a gradient impulse response function (GIRF) [@vannesjo_gradient_2013], is particularly needed for high-quality MRI reconstruction. 

Furthermore, MRI acquisition can be further accelerated by acquiring a larger proportion of k-space per excitation in which each readout duration is prolonged, thus is more susceptible to the impact of static off-resonance (or field inhomogeneity, B<sub>0</sub>). For non-Cartesian imaging, the resulting artifacts can be accounted for by incorporating off-resonance maps into the image reconstruction [@sutton_fast_2003].

This software package, `GIRFReco.jl`, provides an open-source, single ecosystem implementation (in Julia [@bezanson_julia_2017]) of a state-of-the-art image reconstruction for spiral MRI [@wilm_higher_2011;@wilm_diffusion_2015] with a reutilization of core reconstruction implementations in `MRIReco.jl`, a comprehensive open-source image reconstruction toolbox also in Julia. The reconstruction procedure uses a signal model with additional terms represent system imperfections and off-resonance in combination with parallel imaging acceleration and the respective iterative non-Cartesian image reconstruction algorithms (e.g., CG-SENSE [@pruessmann_advances_2001]). In particular, to enable robust and fast MRI with spiral gradient waveforms, `GIRFReco.jl` is designed as an end-to-end signal processing pipeline from raw open-standard MR data ([ISMR]MRD [@inati_ismrm_2017]) to final reconstructed images (NIfTI format commonly used in processing packages for neuroimaging). It integrates system characterization via GIRF for accurate representation of the encoding fields, calibration data (coil sensitivity and static off-resonance maps) and non-Cartesian iterative parallel imaging reconstruction [@vannesjo_image_2016].

# Statement of Need

Currently available open-source implementations for the correction of system imperfections and static off-resonance in MRI are often implemented within the framework of mature image reconstruction suites such as MRecon (https://www.gyrotools.com/gt/index.php/products/reconframe), Gadgetron [@hansen_gadgetron_2013] or BART [@blumenthal_mrireconbart_2022].

However, the aforementioned complexity of the image reconstruction task for spiral MRI currently necessitates the integration of multiple software tools developed in multiple programming languages (C for BART, C++ for Gadgetron, MATLAB for MRecon, etc.) to establish a performant and comprehensive image reconstruction workflow [@veldmann_opensource_2022] that maximize the advantages of each language while minimizing their drawbacks. Extending such pipeline of image reconstruction requires cross-language expertise and knowledge, adding significant overhead and complexity to development. This presents a significant barrier to efficient and reproducible image reconstruction and limits software accessibility and sustainability, especially for non-technical users such as radiologists and other physicians.

The programming language Julia [@bezanson_julia_2017] provides a solution to this multiple-language problem by using a high-level interface to low-level compiled code, i.e., enabling fast prototyping with limited resources in an academic setting, while delivering performant execution times without code rewriting, all within a single developing environment.

In this work, we developed the `GIRFReco.jl`, an image reconstruction pipeline(initial version introduced at the annual meeting of ISMRM 2022 [@jaffray_open-source_2022]) based on the established `MRIReco.jl` package, which implements an end-to-end, self-contained processing and image reconstruction pipeline for spiral MR data completely in Julia, incorporating the different model and correction components [@sutton_fast_2003;@wilm_higher_2011;@wilm_diffusion_2015;@pruessmann_advances_2001;@vannesjo_image_2016] to achieve a high-quality spiral MRI reconstruction. In detail, this reconstruction pipeline consists of several major steps: (1) ESPIRiT coil sensitivity map estimation [@uecker_espirit-eigenvalue_2014]; (2) Robust off-resonance (B<sub>0</sub>) map estimation [@funai_regularized_2008]; (3) GIRF Correction to the theoretical non-Cartesian k-space trajectory [@vannesjo_gradient_2013;@vannesjo_image_2016]; (4) Iterative non-Cartesian MRI reconstruction with off-resonance correction [@knopp_iterative_2009]. Considering software reusabulity and sustainability, (1) and (4) of the abovementioned steps are handeled by `MRIReco.jl`, a comprehensive lower-level open-source image reconstruction toolbox also in Julia. In step (3), the GIRF correction was implemented through our another package [`MRIGradients.jl`](https://github.com/MRI-gradient) [@jaffray_open-source_2022], which is a major refactoring of a MATLAB implementation of the original authors.

<!-- 
In the spirit of software re-usabulity and sustainbility, __GIRFReco.jl__'s core image reconstruction routines are handled by __MRIReco.jl__ [15], a comprehensive open-source image reconstruction toolbox also in Julia. In particular, __GIRFReco.jl__ utilizes __MRIReco.jl__'s iterative non-Cartesian image reconstruction [7], ESPIRiT estimation of coil sensitivity maps [16], and off-resonance correction [17], as well as open-source handling of MR raw data ([ISMR]MRD, [9]). In the course of this work, we implemented robust off-resonance map estimation [18] as a submodule (__MRIFieldmaps.jl__) under __MRIReco.jl__.

With __GIRFReco.jl__, we provide a self-contained, open-source spiral MRI reconstruction pipeline built on __MRIReco.jl__ and include submodules which correct trajectory imperfection using a Julia implementation of the GIRF correction [1,10] (partially refactored from a MATLAB implementation of the original authors https://github.com/MRI-gradient), output corrected MR raw data in ISMRMRD format [9] for the spiral acquisitions, and perform necessary intermediate processing steps. __GIRFReco.jl__ incorporates flexible parameterization of arbitrary input k-space data and supports multi-shot acquisition. It is compatible with GIRF measurements acquired using phantom scans [19,20,21], and thus can provide model-based corrections without need of external hardware [20]. The pipeline executes performant, compiled code, a feature of the Julia language.
-->

# Functionality

## Required Inputs

`GIRFReco.jl` requires raw MR (k-space) data (in [ISMR]MRD format [@inati_ismrm_2017]) of the following scans as input:

1. Multi-echo Gradient-echo spin-warp (Cartesian) scan
    - at least two echo times (e.g., **TODO: xx and yy (water fat in-phase at 3T)**)
2. Spiral scan
    - single or multi-interleave

At the moment, the slice geometries (thickness, field-of-view, and direction) of the Cartesian and spiral scans needs to be matched with each other.

## Overview of Components

The following components are utilized within the spiral reconstruction pipeline of `GIRFReco.jl` (Fig. 1), and separated into different (sub-)packages.

1. Core iterative image reconstruction, using Julia package `MRIReco.jl`
    - CG-SENSE [@pruessmann_advances_2001] algorithm
    - ESPIRiT [@uecker_espirit-eigenvalue_2014] for sensitivity maps
2. Model-based correction components
    - Smoothed B<sub>0</sub> map estimation, custom implementation of [@funai_regularized_2008] in `MRIFieldMaps.jl`
    - Static B<sub>0</sub> map correction, accelerated by time-segmented implementation in `MRIReco.jl` [@knopp_mrirecojl_2021]
    - Gradient impulse response function (GIRF) [@vannesjo_gradient_2013]
        - Measured on a phantom [@graedel_comparison_2017]
        - Open-source computation [@wu_mr_2022]
        - Prediction implemented in our customized package `MRIGradients.jl` [@jaffray_open-source_2022]

![Figure 1](paper/GIRFReco_Components.png?raw=true "GIRFReco.jl Components")

## Detailed Processing Pipeline

`GIRFReco.jl` executes the required steps for spiral diffusion reconstruction (shown in Figure 2) in the following order:

1. Conversion of proprietary format, vendor-specific raw image data to an open-source raw data format ([ISMR]MRD, [@inati_ismrm_2017])
2. Synchronization of the data and the k-space trajectory
3. Model-based correction of the k-space sampling points and data using the gradient impulse response function (GIRF [@vannesjo_gradient_2013], `MRIGradients.jl`)
4. Coil sensitivity map estimation (ESPIRiT, [@uecker_espirit-eigenvalue_2014])
5. Off-resonance (B<sub>0</sub>) map estimation (`MRIFieldmaps.jl`, [@funai_regularized_2008])
6. Non-Cartesian, iterative parallel image reconstruction with off-resonance correction (`MRIReco.jl`, [@knopp_mrirecojl_2021])

![Figure 2](paper/GIRFReco_Pipeline.png?raw=true "Template Reconstruction Pipeline")

## Quality of Life Features

In addition to providing an end-to-end reconstruction workflow, `GIRFReco.jl` provides intermediate methods for plotting images and calibration data using PlotlyJS, as well as the capability to export intermediate reconstruction results (such as calculated coil sensitivity maps and B<sub>0</sub> maps).

# Getting Started

Up-to-date information about how to install `GIRFReco.jl`, run example reconstructions and apply it to your own data can be found in the README.md provided in the GitHub repository. Further technical documentation about the API is provided at **TODO: link to documenter results**, automatically generated by [`Documenter.jl`](https://github.com/JuliaDocs/Documenter.jl).

# Conclusion and Outlook

The presented package, `GIRFReco.jl`, is an open-source end-to-end toolbox for spiral MRI reconstruction. It is developed in Julia, allowing users to obtain final images directly from raw MR data acquired by spiral k-space trajectories. Following best practices of software sustainability and accessibility, we re-use code of the established MR image reconstruction package `MRIReco.jl` in our package, while extending its capability to handle the more complex use case of multiple model-based corrections, necessary for high-quality spiral MRI.

Without loss of generality, besides spiral imaging, `GIRFReco.jl` can be easily extended to a workflow of model-based MR image reconstruction for data acquired under arbitrary non-Cartesian k-space trajectories. With the established functions for GIRF and off-resonance corrections, k-space raw data acquired by other non-Cartesian trajectories (such as radial) will also benefit from high-quality reconstruction by `GIRFReco.jl`.

# References

<!--
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
20. M. Stich et al., "Field camera versus phantom-based measurement of the gradient system transfer function (GSTF) with dwell time compensation,"Magn. Reson. Imaging, vol. 71, pp. 125-131, Sep. 2020, do: 10.1016/j.mri.2020.06.005.
21. R. K. Robison, Z. Li, D. Wang, M. B. Ooi, and J. G. Pipe, "Correction of B<sub>0</sub> eddy current effects in spiral MRI," Magn. Reson. Med., vol. 81, no. 4, pp.2501-2513, 2019, doi:10.1002/mrm.27583.
22. Z. Wu, A. Jaffray, J. Vannesjo, K. Uludag, and L. Kasper, “MR System Stability and Quality Control using Gradient Impulse Response Functions (GIRF),” in Proc. Intl. Soc. Mag. Reson. Med. 30, London, UK, p. 0641. [Online]. Available: https://index.mirasmart.com/ISMRM2022/PDFfiles/0641.html

-->
