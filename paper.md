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
    orcid: 0000-0002-9571-1838
    equal-contrib: True
    corresponding: True
    affiliation: 1
  - name: Zhe Wu
    orcid: 0000-0002-2079-5977
    equal-contrib: True # (This is how you can denote equal contributions between multiple authors)
    corresponding: True # (This is how to denote the corresponding author)
    affiliation: 2
  - name: S. Johanna Vannesjo
    orcid: 0000-0000-0000-0000
    equal-contrib: False # (This is how you can denote equal contributions between multiple authors)
    corresponding: False # (This is how to denote the corresponding author)
    affiliation: 3
  - name: Kamil Uludag
    orcid: 0000-0000-0000-0000
    equal-contrib: False # (This is how you can denote equal contributions between multiple authors)
    corresponding: False # (This is how to denote the corresponding author)
    affiliation: "2, 4" # (Multiple affiliations must be quoted)
  - name: Lars Kasper
    orcid: 0000-0000-0000-0000
    equal-contrib: False # (This is how you can denote equal contributions between multiple authors)
    corresponding: False # (This is how to denote the corresponding author)
    affiliation: "2, 4" # (Multiple affiliations must be quoted)
affiliations:
 - name: MRI Research Centre, University of British Columbia, Vancouver, Canada
   index: 1
 - name: Techna Institute, University Health Network, Ontario, Canada
   index: 2
 - name: Department of Physics, Norwegian University of Science and Technology, Trondheim, Norway
   index: 3
 - name: Department of Medical Biophysics, University of Toronto, Canada
   index: 4
date: 31 January 2023
bibliography: paper/paper.bib

# Optional fields if submitting to a AAS journal too, see this blog post:
# https://blog.joss.theoj.org/2018/12/a-new-collaboration-with-aas-publishing
# aas-doi: 10.3847/xxxxx <- update this with the DOI from AAS once you know it.
# aas-journal: Astrophysical Journal <- The name of the AAS journal.
---

# Summary

Magnetic Resonance Imaging (MRI) acquires data in the frequency domain (k-space), with the sampling pattern traversed by a path known as the k-space trajectory. Traditional Cartesian MRI is inefficient due to intrinsically low signal-to-noise ratio (SNR) per unit time, a consequence of sampling on a discrete grid with short readout durations. On the other hand, non-Cartesian (i.e., non-rectilinear) trajectories with long readouts, such as spirals, can offer considerable improvement in SNR per unit time at the cost of reconstruction complexity [@lee_signalnoise_2021].

The actual k-space trajectory applied during the MRI experiment can differ from the nominal trajectory due to hardware imperfections, resulting in image artifacts such as ghosting, blurring or geometric distortion. This problem is exacerbated in many non-Cartesian trajectories, including spiral readouts, because these fast imaging protocols place high demands on the gradient hardware of the MRI system [@block_spiral_2005].

To combat this, accurate characterization of the system hardware for k-space trajectory correction is necessary, for example via a gradient impulse response function (GIRF) [@vannesjo_gradient_2013]. 

The encoding of the object into the frequency domain is susceptible to static off-resonance (or field inhomogeneity, B<sub>0</sub>), the impact of which is magnified at longer readout durations such as those arising in non-Cartesian acquisitions. However, the resulting artifacts can be accounted for by incorporating off-resonance maps into the image reconstruction [@sutton_fast_2003].

This software package, `GIRFReco.jl`, provides an open-source, single ecosystem implementation (in Julia [@bezanson_julia_2017]) of state-of-the-art spiral image reconstruction [@wilm_higher_2011;@wilm_diffusion_2015]. The core reconstruction routines rely upon `MRIReco.jl`, a comprehensive open-source image reconstruction toolbox also written in Julia. The reconstruction procedure uses an expanded signal model with additional terms which represent system imperfections and off-resonance, in combination with parallel imaging acceleration and iterative image reconstruction algorithms compatible with arbitrary sampling patterns (e.g., CG-SENSE [@pruessmann_advances_2001]). To enable robust, accessible and fast MRI with spiral gradient waveforms, `GIRFReco.jl` is designed as an end-to-end signal processing pipeline, from open-standard raw MR data ([ISMR]MRD [@inati_ismrm_2017]) to final reconstructed images (NIfTI format commonly used in processing packages for neuroimaging). It integrates system characterization information via GIRF correction for accurate representation of the encoding fields, relevant calibration data (coil sensitivity and static off-resonance maps) and non-Cartesian iterative parallel imaging reconstruction [@vannesjo_image_2016].

# Statement of Need

Currently available open-source implementations for the correction of system imperfections and static off-resonance in MRI are often implemented within the framework of mature image reconstruction suites such as MRecon (https://www.gyrotools.com/gt/index.php/products/reconframe), Gadgetron [@hansen_gadgetron_2013] or BART [@blumenthal_mrireconbart_2022].

However, the aforementioned complexity of the image reconstruction task for spiral MRI currently necessitates the integration of multiple software tools developed in multiple programming languages (C for BART, C++ for Gadgetron, MATLAB for MRecon, etc.) to establish a performant and comprehensive image reconstruction workflow [@veldmann_opensource_2022] that maximize the advantages of each language while minimizing their drawbacks. Extending such a pipeline of image reconstruction requires cross-language expertise and knowledge, adding significant overhead and complexity to development. This presents a significant barrier to efficient and reproducible image reconstruction and limits software accessibility and sustainability, especially for non-technical users such as radiologists and other physicians.

The programming language Julia [@bezanson_julia_2017] provides a solution to this multiple-language problem by using a high-level interface to low-level compiled code, i.e., enabling fast prototyping with limited resources in an academic setting, while delivering performant execution times without code rewriting, all within a single developing environment.

In this work, we developed the `GIRFReco.jl`, an image reconstruction pipeline(initial version introduced at the annual meeting of ISMRM 2022 [@jaffray_open-source_2022]) based on the established `MRIReco.jl` package, which implements an end-to-end, self-contained processing and image reconstruction pipeline for spiral MR data completely in Julia. It incorporates model-based corrections [@sutton_fast_2003;@wilm_higher_2011;@wilm_diffusion_2015;@pruessmann_advances_2001;@vannesjo_image_2016] to achieve high-quality spiral MRI reconstructions. In detail, this reconstruction pipeline consists of several major steps: (1) ESPIRiT coil sensitivity map estimation [@uecker_espirit-eigenvalue_2014]; (2) Robust off-resonance (B<sub>0</sub>) map estimation [@funai_regularized_2008;@Lin_Fessler_2020]; (3) GIRF Correction to the theoretical non-Cartesian k-space trajectory [@vannesjo_gradient_2013;@vannesjo_image_2016]; (4) Iterative non-Cartesian MRI reconstruction with off-resonance correction [@knopp_iterative_2009]. Considering software reusabulity and sustainability, (1) and (4) of the abovementioned steps are handeled by `MRIReco.jl`, a comprehensive lower-level open-source image reconstruction toolbox also in Julia. In step (3), the GIRF correction was implemented through the package `MRIGradients.jl` (https://github.com/MRI-gradient) [@jaffray_open-source_2022], which is a major refactoring of a MATLAB implementation of the original authors.

# Functionality

## Required Inputs

`GIRFReco.jl` requires raw MR (k-space) data (in [ISMR]MRD format [@inati_ismrm_2017]) of the following scans as input:

1. Multi-echo Gradient-echo spin-warp (Cartesian) scan
    - at least two echo times (e.g. 4.92 ms and 7.38 ms at 3T)
2. Spiral scan
    - single or multi-interleave

At the moment, the slice geometry (thickness, field-of-view, and direction) of the Cartesian and spiral scans must be congruent.

## Overview of Components

The following components are utilized within the spiral reconstruction pipeline of `GIRFReco.jl` (Fig. 1), and called from their respective packages.

1. Core iterative image reconstruction, using Julia package `MRIReco.jl`
    - CG-SENSE [@pruessmann_advances_2001] algorithm
    - ESPIRiT [@uecker_espirit-eigenvalue_2014] for sensitivity maps
2. Model-based correction components
    - Smoothed B<sub>0</sub> map estimation, custom implementation of [@funai_regularized_2008;@Lin_Fessler_2020] in `MRIFieldMaps.jl`
    - Static B<sub>0</sub> map correction, accelerated by time-segmented implementation in `MRIReco.jl` [@knopp_mrirecojl_2021]
    - Gradient impulse response function (GIRF) [@vannesjo_gradient_2013]
        - Measured on a phantom [@graedel_comparison_2017]
        - Open-source computation [@wu_mr_2022]
        - Prediction implemented in our customized package `MRIGradients.jl` [@jaffray_open-source_2022]

![Figure 1](paper/GIRFReco_Components.png?raw=true "GIRFReco.jl Components")

## Detailed Processing Pipeline

`GIRFReco.jl` executes the steps required for spiral diffusion reconstruction (shown in Figure 2) in the following order:

1. Conversion of proprietary format, vendor-specific raw image data to an open-source raw data format ([ISMR]MRD, [@inati_ismrm_2017])
2. Synchronization of the data and the k-space trajectory onto a common timebase
3. Model-based correction of the k-space sampling points and data using the gradient impulse response function (GIRF [@vannesjo_gradient_2013], `MRIGradients.jl`)
4. Coil sensitivity map estimation (ESPIRiT, [@uecker_espirit-eigenvalue_2014])
5. Off-resonance (B<sub>0</sub>) map estimation (`MRIFieldmaps.jl`, [@funai_regularized_2008;@Lin_Fessler_2020])
6. Non-Cartesian, iterative parallel image reconstruction with off-resonance correction (`MRIReco.jl`, [@knopp_mrirecojl_2021])

![Figure 2](paper/GIRFReco_Pipeline.png?raw=true "Template Reconstruction Pipeline")

## Quality of Life Features

In addition to providing an end-to-end reconstruction workflow, `GIRFReco.jl` provides intermediate methods for plotting images and calibration data using PlotlyJS, as well as the capability to export intermediate reconstruction results such as calculated coil sensitivity maps and B<sub>0</sub> maps to NIfTI format.

# Getting Started

Up-to-date information about how to install `GIRFReco.jl`, run example reconstructions and apply it to your own data can be found in the README.md provided in the GitHub repository. Further technical documentation about the API is provided at https://brain-to.github.io/GIRFReco, automatically generated by [`Documenter.jl`](https://github.com/JuliaDocs/Documenter.jl).

# Conclusion and Outlook

The presented package, `GIRFReco.jl`, is an open-source end-to-end toolbox for spiral MRI reconstruction. It is developed in Julia, and allows users to obtain final images directly from raw MR data acquired by spiral k-space trajectories. Following best practices of software sustainability and accessibility, we rely on the established MR image reconstruction package `MRIReco.jl` in our package, while facilitating the handling of the model-based corrections necessary for high-quality spiral MRI reconstruction.

`GIRFReco.jl` can be extended to handle arbitrary k-space trajectories, and can act as a self-contained template for generalized image reconstruction from raw scan and calibration data to interpretable images. 

# References
