# An open-source pipeline for spiral diffusion MRI Reconstruction in Julia 

## 4 Sentence Summary:

Image reconstruction in magnetic resonance imaging (MRI) requires solving a large linear system which may be ill-posed. Efficient solution to such a system is commonly achieved through the assumption that data points follow a rectilinear grid such that the fast Fourier transform can be utilized to reduce the complexity of the computations. However, it is in many ways advantageous to acquire data in a pattern which deviates significantly from a uniformly dense grid of points. Such a sampling pattern requires additional steps including gridding and re-gridding, as well as compensation for hardware imperfection and physical changes within the imaged subject, including motion and magnetic field effects including eddy currents. 

## Motivation:

Image reconstruction from raw frequency-domain data obtained during an MRI acquisition requires accurate characterization of the MR system. This is especially important in the case of imaging protocols which place high demands on the hardware of the MRI system, or those protocols which are inherently susceptible to external perturbations. Acquisition of diffusion-weighted images with non-Cartesian imaging trajectories requires the application of strong magnetic field gradients for diffusion sensitization, as well as accurate knowledge of the applied gradient sequence to maintain geometric fidelity of the reconstructed image. 

Recent works in MR image reconstruction have provided the research community with a number of open-source tools and algorithms for image reconstruction and correction of imaging artefacts (CITE sigpy, mrireco.jl, bart, b0 correction, girf, etc…). 

## Body:

Building on existing Julia libraries for MR image reconstruction, we developed a comprehensive, open-source pipeline for non-Cartesian image reconstruction including B0 inhomogeneity correction, gradient system characterization and eddy current correction. We demonstrated that spiral images can be reconstructed using this pipeline.

Effective strategies to correct these imperfections have been developed and were
shown to be necessary for accurate spiral imaging reconstruction.2-7

The development of open-source, extensible MR image reconstruction frameworks (MRIReco.jl) capable of handling
arbitrary trajectories and expanded signal models have made efficient and generalized non-Cartesian image reconstruction accessible. 2,8 In the
present work, we demonstrate a comprehensive, open-source pipeline for the reconstruction of images acquired using spiral k-space trajectories.
The pipeline is built around MRIRecojl and ISMRMRD, and incorporates Bo correction, correction for trajectory imperfections using the Gradient
Impulse Response Function (GIRF), and eddy current compensation prior to reconstruction.

The end-to-end reconstruction pipeline from Siemens raw k-space data to spiral images is depicted in Figure 1, highlighting the inclusion of
established open-source 3rd party tools (blue and yellow), as well as the custom modules developed in this work (purple).

Despite these small imperfections, the reconstruction based on the nominal spiral trajectory without Bo correction exhibits severe blurring artifacts, image intensity modulations and geometric inaccuracies (Fig. 4B, 4F). Reconstruction with Bo correction significantly reduces blurring (Fig. 4C, 4G). Reconstruction with the corrected spiral trajectory obtained from GIF prediction reduces geometric inaccuracy due to gradient delays (Fig. 4D, 4H). Finally, incorporation of GIRF-predicted k0 correction (Fig. 3C) leads to virtually artifäct-free spiral image reconstructions with geometric congruency between GRE and spiral scans. (Fig. 4E, 41).

The versatility of the pipeline is shown in Fig. 5 for the reconstruction of b=1000 diffusion-weighted images. The corrections applied are effective
even in the presence of diffusion-weighting gradients, with residual imperfections presumably due to higher-order effects. Reconstruction time with Bo, GIF and k correction was 19 seconds per slice when run on an Apple M1 MacBook Air.

We demonstrated fast, high-quality spiral diffusion imaging enabled by a modular, open-source reconstruction pipeline written in Julia. The
incorporation of GIF and Bo correction was accommodated within Julia and was necessary for accurate spiral reconstruction, as previously
shown 4,11-13 The computational efficiency of Julia and MRIReco.jlenable reasonably fast image reconstruction on modest hardware. In combination
with the open and streaming-capable ISMRMRD format, real-time robust spiral image reconstruction can be realized.° No additional hardware was
needed to achieve high reconstruction quality, and the ISMRMRD file interface makes the framework vendor- and system-agnostic. The modularity
of the presented pipeline lends itself to further development of the reconstruction framework to address remaining imperfections (Fig. 5),2,3
To make such corrections easily accessible to the broader MR imaging community and add to the growing list of MR reconstruction tools available in
Julia,& this reconstruction pipeline will be available on GitHub and as a package within the Julia ecosystem (GIRFReco.jl).
