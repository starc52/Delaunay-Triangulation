# Lawson edge swap implementation 

## Core files:
```dcel/dcel.py```: the file that contains my DCEL implementation. Implementation originally done in the github repo: [dcel](https://github.com/anglyan/dcel) which I adapted to my use-case.
Note that the original implementation is not maintained and has many bugs. 

```ui-tkinter.py```: the file that contains my core ui and approach to the problem. 

## Other files:
```dcel/iodata.py```,```dcel/pyeps.py```, ```dcel/xygrapy.py```: remnants of the original DCEL implementation from the github repo

```dcel.py```: testing the DCEL data structure.

```lawson-edge-swap.py```: testing the edge-swap algorithm on toy, pre-triangluated example ```ex1.txt```

## Usage:
* Install ```numpy```, ```scipy``` and ```tkinter```.
* Open a terminal and go to the directory that contains ```ui-tkinter.py```
* Run ```python ui-tkinter.py```
* Click as many points on the white window as you want. 
* Next press the ```p``` button on the keyboard for partitioning the convex-hull into two x-monotone polygons.
* Next press the ```r``` button on the keyboard to randomly triangulate the set of points. Essentially it triangulates the two x-monotone polygon and merges them.
* Next press the ```d``` button on the keyboard to do the lawson-edge swap algorithm for delaunay triangluation. 

## Comments
Except for the convex hull finding, all implementations were done from scratch. For convex hull, I use ```scipy.spatial.ConvexHull```.
To the best of my knowledge, this algorithm is bug free.

The explanation to my approach is given in the slides and the accompanying video presentation. Further the code is documented for reader's ease of understanding.

Video Presentation: ```cgeom-presentation.mp4```

Slides: ```slide-deck.pdf```
