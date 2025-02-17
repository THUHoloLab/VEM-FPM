<script
  defer
  src="https://cdn.jsdelivr.net/npm/img-comparison-slider@8/dist/index.js"
></script>
<link
  rel="stylesheet"
  href="https://cdn.jsdelivr.net/npm/img-comparison-slider@8/dist/styles.css"
/>

# VEM-FPM

This is the official implementation of variational Expectation Maximization Fourier ptychographic microscopy (VEM-FPM), a novel reconstruction routine for Fourier ptychographic microscopy (FPM). Here, we propose a latent-wavefront physical model for solving Fourier ptychography by introducing a new latent wavefront at the front surface of the image sensor. The inverse problem of FPM is formulated under the framework of variational expectation maximization (VEM). The VEM-FPM alternates between solving a non-convex optimization problem to estimate the latent wavefront in the spatial domain and solving a convex problem to update the sample wavefront and pupil function in the Fourier domain. The VEM-FPM approach enables a stitching-free, full-field reconstruction for Fourier ptychography over a 5.3 mm × 5.3 mm field of view, using a $2.5\times$ objective with a numerical aperture (NA) of 0.08. The synthetic aperture achieves a resolution equivalent to 0.53 NA @ 532 nm wavelength. The execution speed of VEM-FPM is twice as fast as that of state-of-the-art feature-domain [methods](https://opg.optica.org/optica/fulltext.cfm?uri=optica-11-5-634&id=549881) while maintaining comparable reconstruction quality.


<style>
  .slider-example-opacity-and-size {
    --default-handle-width: 100px;
  }

  .slider-example-opacity-and-size:focus {
    --default-handle-opacity: 0;
  }
</style>
<img-comparison-slider class="slider-example-opacity-and-size">
  <img slot="first" width="100%" src="resources/before.png">
  <img slot="second" width="100%" src="resources/after.png">
</img-comparison-slider>
