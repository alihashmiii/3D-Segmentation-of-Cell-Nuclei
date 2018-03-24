# 3D-Segmentation-of-Cell-Nuclei-in-Clusters-
***Segment Cell Nuclei using 3D-Watershed and Quantile Envelopes***

Authors: Jan Liphardt (Stanford) and Ali Hashmi with help from Anton-Antonov for Quantile-Envelopes

Date: 2015

The Mathematica Package `QuantileNucleiSegmentation.m` can be used to segment cell nuclei in spheroid bodies. The algorithm relies on an elaborate 3D watershed segmentation scheme coupled with Quantile Envelopes to bound nuclei.

*Note:* the package is not yet fully optimized for speed. Although I know of a few ways (`SparseArray implementations at a few places` etc..) to speed up the code, I do not have sufficient time at present to implement those schemes.

![Alt-text](https://user-images.githubusercontent.com/10793580/34075250-5417e288-e2c1-11e7-93d7-e2cb29103fd0.png)
