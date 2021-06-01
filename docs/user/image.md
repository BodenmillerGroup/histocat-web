# Image

Image view is used to visualize channels selected in the Channels view as well as segmentation masks.
Zoom/panning functionality is supported.

![Image view](../assets/image-view.png)

Raster images can be exported in TIFF, OME-TIFF or PNG formats by clicking Export button:  

![Image export](../assets/image-export.png)

One useful feature is to get channels statistics for a selected region.
To activate this behavior, please switch **Region** toggle button and then select some area on the raster image using mouse by holding _SHIFT_ key or by selecting `Lasso` tool ![Lasso tool](../assets/lasso-tool-button.png).
As soon as you release left mouse button, you will get general statistics (min/max/mean intensities) for the selected area for all available channels.
This information is displayed in Region panel (see related section in documentation).
