# GIRFReco.jl: An open-source pipeline for spiral MRI Reconstruction in Julia

## Motivation

High-quality image reconstruction from frequency-domain data obtained during an MR image acquisition requires accurate characterization of the system hardware, for example via a gradient impulse response function (GIRF) [1]. If unaccounted for, system imperfections induce image artifacts such as ghosting, blurring or geometric distortion. This is important in the case of fast, SNR-efficient imaging protocols, especially for those employing spiral readout gradient waveforms which place high demands on the hardware of the MRI system [2]. 

Thus, to enable robust and fast MRI with spiral gradient waveforms, an end-to-end signal processing pipeline is required that integrates system characterization via GIRF for accurate representation of the encoding fields with calibration data (coil sensitivity and static off-resonance maps) and non-Cartesian iterative parallel imaging reconstruction.

## Statement of Need

Currently implemented toolboxes for the correction of system imperfection depend upon existing, high-quality, and mature image reconstruction suites, such as MRecon, Gadgetron or BART. However, comprehensive image reconstruction of spiral MR data currently requires simultaneous use of tools from different suites, presenting a barrier to efficient and reproducible image reconstruction. Extending and developing frameworks which are built in multiple languages requires simultaneous expertise across languages and environments, adding significant overhead and complexity to development. Julia provides a means to solve this multiple-language problem by providing a high-level interface to low-level compiled code. 

In the present landscape of open-source image reconstruction toolboxes, MRIReco.jl facilitates execution of the described reconstruction pipeline by providing non-Cartesian image reconstruction [2020 Knopp], sensitivity map estimation using ESPiRIT [2013 Uecker], robust off-resonance correction [2009 Knopp], and open-source data handling [2011 Kozerke et al.]. 

With GIRFReco.jl, we provide a self-contained, open-source spiral diffusion reconstruction pipeline built on MRIReco.jl and include tools which correct trajectory imperfection using a Julia implementation of GIRF correction [Vannesjo mrigradients], validate ISMRMRD data acquired using spiral trajectories, and perform necessary intermediate processing steps. GIRFReco.jl is entirely written in the Julia language, and interfaces with the ISMRM Raw Data format. It incorporates flexible parameterization of arbitrary input k-space data and supports multi-shot acquisition. It is compatible with GIRF measurements acquired using phantom scans, and thus can provide model-based corrections without the need of external hardware.

## Functionality

The required steps for spiral diffusion reconstruction (shown in Figure 1) are the following:

1.	Conversion of proprietary format, vendor-specific raw image data to an open-source raw data format (ISMRMD)
2.	Synchronization of the data and the k-space trajectory
3.	Model-based correction of the k-space sampling points and data using the gradient impulse response function (GIRF) 
4.	Coil sensitivity map estimation
5.	Off-resonance (B0) map estimation 
6.	Non-Cartesian parallel image reconstruction with off-resonance correction

## Conclusion

We provide here an example reconstruction pipeline that is built on an established reconstruction package (MRIReco.jl) and extend its capability to provide an end-to-end template for spiral image reconstruction which is open-source from raw data to final image. 


## References:

1. S. J. Vannesjo et al., "Gradient system characterization by impulse response measurements with a dynamic field camera," Magn. Reson. Med., vol. 69,no. 2, pp. 583-93, Feb. 2013, doi: 10.1002/mrm.24263.
2. K. T. Block and J. Frahm, "Spiral imaging: A critical appraisal," J. Magn. Reson. Imaging, vol. 21, no. 6, pp. 657-668, Jun. 2005, doi: 10.1002/jmri.20320
2. B. J. Wilm, C. Barmet, M. Pavan, and K. P. Pressman, "Higher order reconstruction for MRI in the presence of spatiotemporal field perturbations:
Higher Order Reconstruction for MRI," Magn. Reson. Med., vol. 65, no. 6, pp. 1690-1701, Jun. 2011, doi: 10.1002/mrm.22767.
3. B. J. Wilm et al., "Diffusion MRI with concurrent magnetic field monitoring," Magn. Reson. Med., vol. 74, no. 4, pp. 925-933, 2015, doi:
10.1002/mrm.25827.
4. S. J. Vannesjo et al., "Image reconstruction using a gradient impulse response model for trajectory prediction: GIRF-Based Image Reconstruction," Magn.Reson. Med., vol. 76, no. 1, pp. 45-58, Jul. 2016, doi: 10.1002/mrm.25841.
5. A. K. Funai, J. A. Fessler, D. T. B. Yeo, V. T. Olafsson, and D. C. Noll, "Regularized Field Map Estimation in MRI," IEEE Trans. Med. Imaging, vol. 27, no.10, pp. 1484-1494, Oct. 2008, doi: 10.1109/TMI.2008.923956.
6. M. Uecker et al., "ESPIRIT-an eigenvalue approach to autocalibrating parallel MRI: Where SENSE meets GRAPPA," Magn. Reson. Med., vol. 71, no. 3,pp. 990-1001, Mar. 2014, doi: 10.1002/mrm.24751.
7. T. Knopp and M. Grosser, "MRIReco./I: An MRI reconstruction tramework written in Julia," Magn. Reson. Med., vol. 86, no. 3, pp. 1,633-1646, Sep. 2021,doi:10.1002/mrm.28792.
8. S. J. Inati et al., "ISMRM Raw data format: A proposed standard for MRI raw datasets," Magn. Reson. Med., vol. 77, no. 1, pp. 411-421, Jan. 2017, doi:10.1002/mrm.26089.
9. M. Lustig, S.-J. Kim, and J. M. Pauly, "A fast method for designing time-optimal gradient waveforms for arbitrary k-space trajectories," IEEE Trans. Med. Imaging, vol. 27, no. 6, pp. 866-873, Jun. 2008, doi: 10.1109/TM1.2008.922699.
10. N. N. Graedel, S. A. Hurley, S. Clare, K. L. Miller, K. P. Pressmann, and S. J. Vannesjo, "Comparison of gradient impulse response functions measured with a dynamic field camera and a phantom-based technique," Barcelona/ES, 2017, p. 378.
11. M. Stich et al., "Field camera versus phantom-based measurement of the gradient system transfer function (GSTF) with dwell time compensation,"Magn. Reson. Imaging, vol. 71, pp. 125-131, Sep. 2020, do: 10.1016/.mri.2020.06.005.
12. R. K. Robison, Z. Li, D. Wang, M. B. Ooi, and J. G. Pipe, "Correction of B0 eddy current effects in spiral MRI," Magn. Reson. Med., vol. 81, no. 4, pp.2501-2513, 2019, doi:10.1002/mrm.27583.
