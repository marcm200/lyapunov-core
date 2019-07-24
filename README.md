# lyapunov-core
Text-menu based code to compute Lyapunov-value based images by Mario Markus' algorithm
from Spektrum der Wissenschaft 1995

# Purpose of this code

To calculate large Lyapunov diagrams and automatically explore parameters in a structured or randomized fashion.

This readme is organized as:

1. Quick start
2. Background
3. Hard-coded functions to use as f or g
4. Commands of the program
5. Limitations
6. Contact
7. Description of the *.par file entries



## (1) Quick start

Compile the program in a C++ compiler, best with speed and math optimization flags on, and run it via (assuming the file is named lyapunov.exe throughout this document).

`lyapunov.exe <quickstart.txt`

This calculates the images whose parameters are stored in the files named `NNtemplate.par` in a small 600x600 version and saves them as bitmaps.

To see the automatic exploration process, run:

`lyapunov.exe <walk.txt`

The results are stored under the files beginning with `_walk...`


## (2) Background

In 1995, Mario Markus published an article in the german science magazine Spektrum der Wissenschaft explaining the use of the Lyapunov value to generate images using a disturbance parameter r that alters the argument of a function and oscillates between the coordinates of a complex point. That oscillation is controlled by a sequence of letters, e.g. AAABBBAAABBB..., which extends indefinitely and is used during the iterating process of computing the Lyapunov exponent for a pixel, wherein A puts the x coordinate of said pixel (in complex plane coordinate) into r and B the y coordinate.

To create a Ljapunow image, one needs two functions, here called f and g. f plays the role of the trajectory (iterating) function which starts from a value x=x0=0.5 and performs the iteration x1=f(x0), x2=f(x1) and so on. After a predefined number of such steps (skipping steps) the actual computation of the Lyapunov exponent starts. The iteration is continued. In each step the value g(x<sub>n</sub>) is computed and log(|g(x<sub>n</sub>|) is added to a sum and at the end averaged to give the final value. That function g is usually the derivative f' of f (but does not need to be, see below).

The software uses hard-coded function pairs f,g as an object (Function_NNN). The coloring here is achieved by defining one to several intervals of Lyapunov values, e.g. [0.1..0.2] and assigning a left and a right RGB value, linearly interpolating between the boundary values if a computed Lyapunov exponent falls within the interval. There are also two intervals from -infinity to a specific value and from another specific value to +infinity which are covered by a constant color.


## (3) Hard-coded functions to use as f or g

Notation below: r is always the disturbance parameter, x is the iterated variable, b is a real number used to add some variation to the function definition to explore (and not to be confused with the sequence character B).

General functions provide procedures called `èval` which compute just the trajectory or the trajectory and the value of g, depending on the number of arguments given.

There is a meta function object which can be used to construct a sectionally defined function pair which - depending on the value of the iterated variable x - uses one of two predefined simple functions f1 and f2 and their corresponding computing functions g1 and g2. That often gives rise to new behavior and images. The file `17template.par` demonstrates the use of that construction. The section definition can be different for the initial skipped iterations (variables i0min and i0max in the parameter file) and the later computing ones (variable i1min and i1max).

E.g. initial 500 iterations: f(x)=r * x if 0 <= x <= 1, else f(x)=sin(x+r)
<br>computing 1000 iterations: f(x)=r * x if -1 <= x <= 0.5, else f(x)=sin(x+r)

Another generilization of the f,f' combination is allowing any function whatsoever to take on the role of g, not necessarily the derivative of f. Some functions are hard-coded that way. But there is also a meta function object called Meta detached, which takes two simple function objects h1 and h2. It can be defined whether meta.f is h1.f or h1.g and whether meta.g is h2.f or h2.g. That way one can combine every function as a trajectory with every function as a computing one and create new behaviour as in the sectionally defined case. The file `16template.par` demonstrates the use of that.

The following functions are currently supported - internally represented by a small integer:

<table>
  <tr><td>Internal number</td><td>definition</td></tr>
  <tr><td>1</td><td>logistic equation f(x)=r*x+(1-x), g=f'</td></tr>
    <tr><td>2</td><td>f(x)=b*sin^2(x+r), g=f'</td></tr>
  <tr><td>3</td><td>f(x)=b*sin(x+r*cos(x+r)), g=f'</td></tr>
<tr><td>7</td><td>f(x)=b*sin(x+r)*sin(x-r), g=f'</td></tr>
<tr><td>11</td><td>f(x)=b*sin^2(x+r), DETACHED g=r-2rx <> f'</td></tr>
<tr><td>13</td><td>f(x)=b*sin(x+r)+b*sin^2(b*x+r), DETACHED g=sin^2(x+r*b)-r*x <> f'</td></tr>
<tr><td>14</td><td>f(x)=r*sin^2(x-r)+b*sin^3(x+2r), DETACHED g=rx-b*sin^4(rx-b) <> f'</td></tr>
<tr><td>16</td><td>Meta function object for defining a detached function pair f,g</td></tr>
<tr><td>17</td><td>Meta function object for sectionally defining a function pair f,g</td></tr>
<tr><td>18</td><td>f(x)=r*sin(x)*(1-b*sin(x+r)), g=f'</td></tr>
  <tr><td>22</td><td>f(x)=b*atan((x+r)*sin(x+r)), g=f'</td></tr>
</table>

## (4) Commands in the text-based menu

Commands are case-insensitive. White spaces should be avoided.

### Main commands

<table>
  <tr><td>RUN</td><td>Calculates the current image with the currently displayed parameters (image size, sequence, position in the complex plane etc). At the end data is automatically stored under the file name tmpljap (see below, storing).</td></tr>

<tr><td>RUN(a,b)</td><td>Calculates only the rows from [a..b] starting from 0 <= a as the bottom row and a,b < height of image. This can be used to split the calculation process in several parts, storing the already computed raw data, reloading it and continue the computation process another time.</td></tr>

<tr><td>E</td><td>Exits the program</td></tr>
</table>


### Loading and saving

<table>
<tr><td>LOAD(filename)</td><td>Filename should be without extension. The command loads filename.par  which stores in a text-manner the parameters of an image. If present, the calculated Lyapunov exponents are loaded from the identically named *.ljd file.</td></tr>

<tr><td>LOADCOLOR(filename)</td><td>Loads only the color method from a given parameter file and overwrites the current color method in memory. If already computed Ljapunow exponents are present in memory, using the save command below generates the new bitmap.</td></tr>

<tr><td>SAVE(filename)</td><td>Stores the parameter file (filename.par), the exponents (extension .ljd), the 24 bit bitmap (.bmp) and a textual description of the image including color method (.descr).</td></tr>
</table>


### Setting some parameters

The parameters are usually set in the parameter file via a simple text editor. Some values can be altered directly:

<table>
<tr><td>SETITER(iter0,iter1)</td><td>The next calculation now uses iter0 (integer) skipping iterations before the start of actual Lyapunov exponent calculation in then iter1 iterations.</td></tr>

<tr><td>SETSIZE(x,y)</td><td>Sets the image size to x columns (integer) and y rows (rounded towards the nearest smaller value divisible by 4).</td></tr>

<tr><td>SETSEQUENCE(string of As and Bs)</td><td>Currently only the 2D version with symbols A and B is supported. </td></tr>

<tr><td>SETPOSITION(a,b,c,d,e,f)</td><td>Sets the position of the rhomboid to lower left (a,b), lower right (c,d) and upper left (e,f). Values are complex plane coordinates.</td></tr>
</table>


### Rhomboid manipulation

The rhomboid that is used to represent the region of the plane to be computed is defined by three complex points: lower left, lower right and upper left. Those values can  be altered via:

<table>
<tr><td>STRETCH(factor)</td><td>Stretches or shrinks it by the given factor (real number) around the center point, thereby covering more or less ground in the complex plane.</td></tr>

<tr><td>ROTATEDEG(value)</td><td>Rotates the rhomboid around the center point via value (real) degrees (360 being a full turn).</td></tr>

<tr><td>CENTER(pixel x,pixel y)</td><td>Centers the current rhomboid to the complex value represented by the pixel coordinates given (as e.g. determined by displaying the picture in an image viewer).</td></tr>

<tr><td>CROP(lower left pixel x,pixel y,upper right pixel x,pixel y)</td><td>Crops the rhomboid via pixel coordinates as given (see CENTER).</td></tr>
</table>


### Automatic exploration

Example usage can be found in walk.txt which can be used as a skript file (see quick start).

<table>
<tr><td>WALKTILE(n,m)</td><td>Tiles the current rhomboid in n tiles horizontally and m tiles vertically and computes every one of them with the current settings (i.e. iteration number, color, image size etc.) and saves those seperately as `_walktile*` files.</td></tr>

<tr><td>WALKRGB</td><td>Generates random RGB values to be put into the current color method (not changing the interval limits though) and saving the parameters and images. Additionally, if the subdirectory `colorcollection\` is present, for every `*.par` file therein its color method is loaded and applied to the current Lyapunov exponents in memory.<br><b>NOTE: This command does not compute the image anew. It uses already available Lyapunov values in memory.</b></td></tr>

<tr><td>WALKB(c,d,n)</td><td>Since almost every function coded has a parameter b, this is now iterated from [c..d] in equally spaced steps and an image is computed with the current loaded settings. Images, parameters and exponents are saved under `_walkb_*` files.</td></tr>

<tr><td>WALKSEQ(n,length)</td><td>Generates randomly a number of sequences of the given length, calculates the images with the current settings and saves them under `_walkseq_*`.</td></tr>

<tr><td>WALKSECTION</td><td>Only appropriate for the sectionally defined Meta function object. The section parameters i0min, i0max, i1min, i1max are iterated between -1..+1 in a small number of equally spaced steps and any combination. Images are computed with the current settings and saved under `_walksection_*`.</td></tr>

<tr><td>WALKDET(func,der0,der1,b0,b1,n)</td><td>This creates a new function of type Meta detached. Its trajectory function part is set to function number func (integer value, see above list) and uses both func.f and func.g as meta.f, its function part g is iterated between der0 and der1 (integer values from above list, not supported numbers are ignored) and also uses that .f and .g part as meta.g. The b value is iterated from b0 to b1 in n steps. Images are computed with the current settings and saved under `_walkdet_*`.</td></tr>
</table>


## (5) Limitations

- The software comes without any warranty.
- It is designed to compute the images. Manual parameter alterations have to be done on the definition file `*.par` in a text editor and for the pixel coordinates an image viewer.
- Functions are hardcoded except for the Meta-functions which allow for arbitrary combinations of hard-coded functions at the cost of lower speed.
- Currently vectorization is not employed.
- There is no special error handling other than simple error messages.
- The bitmap data type was not thoroughly tested to save viewable images of any arbitrary size, but used for images whose size is quadratic and a power of 2 or some easy values like 600 or 800.


## (6) Contact

Please direct comments to: marcm200@freenet.de

Marc Meidlinger, July 2019

Forum: https://fractalforums.org/image-threads/25/lyapunov-diagrams/2273


## (7) Description of the parameter file entries

The `*.par` file consists of one entry per line, case-insensitive, but white-spaces should not be used outside comments (lines starting with #). The general structure is one main tag per line followed by its parameters, which in themselves can consist of several sub tags and the actual parameter values.

<table>
<tr><td>Entry</td><td>Type</td><tdUsage</td></tr>
<tr><td>FUNKTION</td><td>main tag</td><td>description of function to be computed starts here</td><tr>
<tr><td>ID</td><td>sub tag</td></tr>
<tr><td>14</td><td>parameter</td><td>positive integer value that internally represents a given function</td></tr>
<tr><td>B</td><td>sub tag</td><td>this function takes a parameter (not to be confused with the sequence B character)</td></tr>
<tr><td>2.800000e+000</td><td>parameter</td><td>real number</td></tr>
<tr><td>FAERBUNG</td><td>main tag</td><td>color method starts here</td></tr>
<tr><td>ID</td><td>sub tag</td></tr>
<tr><td>2</td><td>parameter</td><td>integer value representing the chosen color method (currently only one implemented)</td></tr>
<tr><td>FARBEL</td><td>sub tag</td><td>if Lyapunov value lies left to every color interval, use the following constant color</td></tr>
<tr><td>255</td><td>parameter</td><td>integer value 0..255 of RGB's red channel</td></tr>
<tr><td>255</td><td>parameter</td><td>integer value 0..255 of RGB's green channel</td></tr>
<tr><td>255</td><td>parameter</td><td>integer value 0..255 of RGB's blue channel</td></tr>
<tr><td>FARBER</td><td>sub tag</td><td>color if Lyapunov value lies right to every color interval starts here</td></tr>
<tr><td>255</td><td>parameter</td><td>as above</td></tr>
<tr><td>255</td><td>parameter</td><td>s above</td></tr>
<tr><td>255</td><td>parameter</td><td>s above</td></tr>
<tr><td>INTANZ</td><td>sub tag</td><td>number of </tr>
<tr><td>2</td><td>parameter</td><td>positive integer denoting the number of compact color intervals</td></tr>
<tr><td>GRENZEL</td><td>sub tag</td><td>left border Lyapunov value for that color interval</td></tr>
<tr><td>-1.916722e+000</td><td>parameter</td><td>real number</td></tr>
<tr><td>GRENZER</td><td>sub tag</td><td>right border</td></tr>
<tr><td>-1.000000e+000</td></tr>
<tr><td>FARBEL</td><td>sub tag</td><td>left border color</td></tr>
<tr><td>0</td></tr>
<tr><td>0</td></tr>
<tr><td>0</td></tr>
<tr><td>FARBER</td><td>sub tag</td><td>right border color</td></tr>
<tr><td>255</td></tr>
<tr><td>255</td></tr>
<tr><td>255</td></tr>
<tr><td>GRENZEL</td><td>sub tag</td><td>next color interval starts here</td></tr>
<tr><td>-1.000000e+000</td><td>parameter</td><td>left border</td></tr>
<tr><td>GRENZER</td><td>as above</td></tr>
<tr><td>0.000000e+000</td></tr>
<tr><td>FARBEL</td></tr>
<tr><td>0</td></tr>
<tr><td>0</td></tr>
<tr><td>0</td></tr>
<tr><td>FARBER</td></tr>
<tr><td>255</td></tr>
<tr><td>255</td></tr>
<tr><td>255</td></tr>
<tr><td>LENX</td><td>main tag</td><td>width of the image</td></tr>
<tr><td>600</td><td>parameter</td><td>positive integer</td></tr>
<tr><td>LENY</td><td>main tag</td><td>height of the image</td></tr>
<tr><td>600</td><td>parameter</td><td>positive integer</td></tr>
<tr><td>ITER0</td><td>main tag</td><td>skipping interation count</td></tr>
<tr><td>250</td><td>parameter</td><td>positive integer</td></tr>
<tr><td>ITER1</td><td>main tag</td><td>computing iteration count</td></tr>
<tr><td>500</td><td>parameter</td><td>positive integer</td></tr>
<tr><td>X0</td><td>main tag</td><td>starting x0 value for iteration</td></tr>
<tr><td>5.000000e-001</td><td>parameter</td><td>real number</td></tr>
<tr><td>SEQUENZ</td><td>main tag</td><td>disturbance sequence starts here</td></tr>
<tr><td>AAAAABBABBBB</td><td>parameter</td><td>character string of A and B in any combination</td></tr>
<tr><td>OL</td><td>main tag</td><td>upper left corner of the rhomboid starts here</td></tr>
<tr><td>-3.200000e+000</td><td>parameter</td><td>real number denotes x coordinate in complex plane</td></tr>
<tr><td>3.200000e+000</td><td>parameter</td><td>real number denotes y coordinate</td></tr>
<tr><td>UL</td><td>main tag</td><td>lower left corner point</td></tr>
<tr><td>-3.200000e+000</td></tr>
<tr><td>-3.200000e+000</td></tr>
<tr><td>UR</td><td>main tag</td><td>lower right corner point</td></tr>
<tr><td>3.200000e+000</td></tr>
<tr><td>-3.200000e+000</td></tr>
  </table>
`

