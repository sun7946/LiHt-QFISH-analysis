# LiHt-QFISH-analysis
1. The file "dot_signal_TimelessM.m" contains the functions used for image processing and signal extraction. 
2. The file "zstack_prepro_looper.m" is a running script for correcting and segmenting the images of the telomere channel. Its input is an .nd2 file containing the z-stack, and its outputs include: ① the maximum intensity projection image; ② the image after removing the camera CMOS, the bias of light source , and the background of the culture medium; ③ the image after removing the internal background of the cells; ④ the mask image of the telomeres.
3. The file "Telemere_fixed_intensity_looper.m" is a running script for extracting the signals of the telomere channel and other channels of interest. Its inputs are: ① an .nd2 file containing the z-stack; ② the mask image of the nuclear channel; ③ the mask image of the telomeres; ④ the output ③ of the "file zstack_prepro_looper.m". Its output is the information of all the cells in each field of view.  
The output contains the following content:  
   x: The x-coordinate of the cell nucleus.  
   y: The y-coordinate of the cell nucleus.  
   area: The area of the cell nucleus.  
   ch1_median: Median intensity of channel 1 in cell nuclear.  
   ch1_ring_median：Median intensity of channel 1 in cell cytoplasmic loop.  
   ch1_sum: Sum intensity of channel 1 in cell nuclear.  
   ch1_mean: Mean intensity of channel 1 in cell nuclear.  
   ch_telomere_mean: Mean intensity of each telomere in cell nuclear.  
   ch_telomere_median: Median intensity of each telomere in cell nuclear.  
   ch_telomere_sum: Sum intensity of each telomere in cell nuclear.  
   ch_telomere_dot_Centroid: The xy-coordinate of each telomere.  
   ch_telomere_dot_area: The area of each telomere.  
   ch_telomere_dot_median_in_median: Median intensity of all telomere's median intensity in a single cell nuclear.  
   ch_telomere_dot_median_in_sum: Median intensity of all telomere's sum intensity in a single cell nuclear.  
   ch_telomere_dot_sum_in_median: Sum intensity of all telomere's median intensity in a single cell nuclear.  
   ch_telomere_dot_sum_in_sum: Sum intensity of all telomere's sum intensity in a single cell nuclear.  
   ch_telomere_dot_number: The number of telomere in a single cell nuclear.
