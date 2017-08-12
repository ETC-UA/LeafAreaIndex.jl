# Gap Fraction

## Binarization

In a first step each pixel of the image is categorized as either sky or vegetation. This is called binarization. 

The binarization methods in this package all require the ability to order image pixels (technically, to be able to compare two pixels and determine if one's values is less than the other's).
Therefore we work on the basis of a single light intensity value for each pixel, as this is not defined unambiguously for RGB values. We recommend using the Image package to select the  blue channel using the `blue` method. If you wish, however, you can easily use any linear combination of the different colour channels.

In general, these binarization methods determine a single global threshold between sky and vegatation. Any pixel with values above this threshold is considered sky.

The standard method is from [Ridler Calvard, 1978](http://ieeexplore.ieee.org/stamp/stamp.jsp?arnumber=4310039), which is similar to [k-means clustering](http://en.wikipedia.org/wiki/K-means_clustering) with two classes in one dimension.

The second method is the [Minimum method](ftp://krill.antcrc.utas.edu.au/pub/el2/Papers/Thresholding/Glasbey,%201993.pdf) that looks for the minimum between two peaks in the histogram.

Finally, the Edge Detection method by [Nobis & Hunziker, 2004](http://www.slf.ch/info/mitarbeitende/nobis/veroeffentlichungen_EN/publications/2005_NOBIS_HUNZIKER.pdf) aims to maximize the contrast at edges between sky and vegetation.

A comparison of different methods can be found in [Glatthorn & Beckschäfer, 2014](http://journals.plos.org/plosone/article?id=10.1371/journal.pone.0111924).

You can access the different methods with (where `polimg` is a `PolarImage` instance):

* `thresh_RC = threshold(polimg)` or `thresh_RC = RidlerCalvard(polimg)`
* `thresh_min = minimum_threshold(polimg)`
* `thresh_edge = edge_threshold(polimg)`

The `RidlerCalvard` and `minimum_threshold` only use pixels within a 90ᵒ view angle, while `edge_threshold` cuts out a rectangle around the visible pixels to maintain order of neighbouring pixels.



## Gap Fraction

Once a threshold ahs been set, you can calculate a gap fraction for a specific part of the image and specialized methods exist for its logarithm and the contact frequency. However, these low level functions are in general not required to be called by the user but are automatically called by higher level functions. See the implementation section for more details.