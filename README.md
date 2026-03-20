# Myofascial-Pain-Project-
Use neck data to quantify short T2* components from fast UTE acquisitions, with the goal of identifying potential biomarkers for MPS.


Location- Source file: Desktop/randika-MRI data/1.stanford_data/niftis_anon(dataset01)/
Location- Work files: Desktop/randika-MRI data/my work/1.stanford_data/1.new_beginning/1.new_beginning.......

In this UTE dataset, no phase data is available, so fat suppression is performed directly using a three-compartment model.
Dixon water–fat separation is not required, as all components are accounted for within the multi-compartment fitting.

Preprocessing & Registration-------------------------------------------------------

A muscle map is required for segmentation.

Neck water images are used as the reference for:

Registration

Reslicing between UTE and water datasets

The resulting transformation (from Dixon/water data) is then applied to generate the muscle segmentation maps.


Component Separation-----------------------------------------------------------------------------

After segmentation, run the Rician-based three-compartment fitting directly on the UTE data to separate:

Short T2* components

Long T2* components

Fat fraction

Spatial Alignment Fix---------------------------------------------------------------------------------

There is a mismatch between output maps and source images (not in the same space).
This is corrected using FSL:

fslcpgeom merged_ute.nii.gz T2str_water_long_map.nii.gz

Quantification---------------------------------------------------------------------------

A separate mask directory contains:

ROI masks (including eroded masks)

Scripts for quantification

These scripts compute:

Short T2* fraction

Long T2* fraction

Fat fraction

Outputs are saved as Excel files for downstream analysis.

Validation---------------------------------------------------------------------------------

Dixon-based water and fat separation outputs were used to compute fat fraction (FF), which was then compared against UTE-derived FF measurements to validate the results.

Notes------------------------------------------------------------------------------------

Ensure all images are properly aligned before quantification.

Eroded masks are recommended to minimize partial volume effects



