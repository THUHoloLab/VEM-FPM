# Latent wavefront model for Fourier ptychography

<div align="center">
<img src="resources/EM_algo.gif" width = "830" alt="" align=center />
</div><br>

This is the official implementation of variational Expectation Maximization Fourier ptychographic microscopy (VEM-FPM), a novel reconstruction routine for Fourier ptychographic microscopy (FPM). Here, we propose a latent-wavefront physical model for solving Fourier ptychography by introducing a new latent wavefront at the front surface of the image sensor. The inverse problem of FPM is formulated under the framework of variational expectation maximization (VEM). The VEM-FPM alternates between solving a non-convex optimization problem to estimate the latent wavefront in the spatial domain and solving a convex problem to update the sample wavefront and pupil function in the Fourier domain. The VEM-FPM approach enables a stitching-free, full-field reconstruction for Fourier ptychography over a 5.3 mm × 5.3 mm field of view, using a $2.5\times$ objective with a numerical aperture (NA) of 0.08. The synthetic aperture achieves a resolution equivalent to 0.53 NA @ 532 nm wavelength. The execution speed of VEM-FPM is twice as fast as that of state-of-the-art feature-domain [methods](https://opg.optica.org/optica/fulltext.cfm?uri=optica-11-5-634&id=549881) while maintaining comparable reconstruction quality.

## How to use
MATLAB codes for VEM-FPM are available in the folder "VEM_FPM". 

Simply run "VEM_FPM_main.m" to begin reconstrucion. The codes were tested on a desktop, running Windows 11 Pro OS, with a CPU of Intel Core i9-12900K, 3.2Hz; a GPU of NVIDIA GeForce RTX 3090.

## Variational Expectation Maximization
The following image shows the working pipeline for VEM-FPM. (a) Optical layout of a Fourier ptychographic microscopy system. (b) Sketch for variation EM algorithm. The VEM alternatively finds the estimation of expectation denoted by the green parabola ($E$-step), and maximizes the likelihood by locating the lateral coordinate of the vertex of the parabola ($M$-step). (c) is the flowchart of VEM-FPM. The latent variable $\mathbf{Z}$ is estimated in the $E$-step through the observed data and predicted data generated by the forward model. In the $M$-step, the sample wavefront $\mathbf{\Phi}$ and the pupil function $\mathbf{P}$ are obtained through blind deconvolution by solving the quadratic loss function with plug-and-play prior for in-loop image denoising. 
<div align="center">
<img src="resources/layout.png" width = "760" alt="" align=center />
</div><br>

## Comparison results
FPM reconstruction using VEM-FPM. (a) full-field reconstruction results and the first 9 images for the red channel. (b) quantitative intensity profile along the white lines in (c) and (d). (c, d) Zoom-in image for the yellow boxes in (a). Scale bar: 500 um for (a); 100 um for (c) and (d). A full-resolution image is available on [[Gigapan]](http://gigapan.com/galleries/13891/gigapans/235892).

<div align="center">
<img src="resources/results.jpg" width = "730" alt="" align=center />
</div><br>

## Benchmarks
The performance of VEM-FPM is compared to state-of-the-art FPM reconstruction methods. Benchmarking can be found at [[FPM benchmark]](https://github.com/ShuheZhang-MUMC/FPM_benchmarks), using simulation data and published data. 
