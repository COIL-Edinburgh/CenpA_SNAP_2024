# Centromere_Quantifications_JPLab

Plugins for quantification of centromeric signal adapted from the CRaQ Macro (1) 

INSTALLATION

Save the CenpA_SNAP_2024.class and/or EndogenousCentromeres_Intensity.class file in the plugins folder of your ImageJ installation. ImageJ should be version 1.5 or later and should have BioFormats installed. The plugin should appear in the Plugins menu.

USAGE

- Select the plugin, if you want to analyse just one (Data) channel use CenpA_SNAP_2024 or two data channes use EndogenousCentromere_Intensity.
- The plugin will ask for an input folder and open all .tif or .ims files in the selected folder
- If the opened file has 4 channels a Maximum Z-Projection is performed. For the first file in the folder a dialog box will
appear asking for the channel number (0-3) for the Data, Reference, DAPI and other Channels. There is the option at this stage to edit the default parameter settings.
- Default parameter settings are:
  - Square Size = 7 pixels
  - Minimum Circularity = 0.95 au
  - Max Feret Diameter = 7.0 pixels
  - Min Centromere Size = 4 pixels
  - Max Centromer Size = 35 pixels
  - No chromatic abberation offsets
- The plugin then for each file:
  - Creates a mask of the DAPI channel to identify cell ROIs
  - Creates a mask of the Reference Channel to identify Centromere ROIs (using parameter settings to filter centromere ROIs)
  - Draws a Square (of size set in parameters) ROI around each Centromere ROI and returns the cell the roi is associated with and the max-min intensity in that region in the data channel (and in the 'other' channel for the EndogenousCentromere_Intensity plugin).

OUTPUT
- An OUTPUT.csv file with the cell and intensity data for each image analysed
- A Log file with the parameter settings and Channel order
- Annotated .tif images showing the ROIs analysed

  (1) Bodor DL, Rodríguez MG, Moreno N, Jansen LET. 2012. Analysis of protein turnover by quantitative SNAP-based Pulse-Chase imaging. Current Protocols in Cell Biology Chapter 8, Unit8.8. doi: 10.1002/0471143030.cb0808s55
