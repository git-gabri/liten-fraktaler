# Liten fraktaler
A fractal rendering software, written in C++, inspired by a much more popular program called [Kalles Fraktaler](https://mathr.co.uk/kf/kf.html), and [FractalC](https://github.com/ahhhh6980/FractalC).

## Compiling on Linux with `gcc`
* Clone this repository
* `cd` into the repository
* `g++ -o liten_fraktaler main.cpp -Wall -O2 -lpthread -lpng`

You should now have an executable called `liten-fraktaler`.  
A compiled (with `gcc 7.5.0`) executable for `x86-64` is also provided in the `bin` folder.

## Usage
This help is also available by launching the program with the `-H` flag, like this `./liten-fraktaler -H`.

Usage: `./liten-fraktaler [OPTIONS]`  
Options:
* **Help**
	* `-H` show this help  
* **Image related flags**
	* `-w SIZE_T` sets width of the output image. The default is 1920.  
	* `-h SIZE_T` sets height of the image. The default is 1080.
	* `-S SIZE_T` sets maximum sector side.  
		To utilize multiple processors to speed up the rendering this program subdivides the image in sectors, whose dimensions are at maximum `256x256` pixels, by default. This flag allows the user to change the maximum side of the sectors.
    * `-o STRING` sets output file name to the content of `STRING`. The default name is `fractal.png`.  
    	Do not include the extension, the output image is always a PNG file.
* **Fractal related flags**
	* `-f INT` sets fractal type. Must be a number in [0, 3], the default is 0.  
		* `0` - Mandelbrot set  
        * `1` - Burning ship  
        * `2` - Mandelbar  
        * `3` - Spade fractal (z^z + c), credits to [ahhhh6980](https://github.com/ahhhh6980)
    * `-r LNG_DBL` sets real part of the offset from the origin
    * `-i LNG_DBL` sets imaginary part of the offset from the origin  
    	The rendering of the fractal is done around a point whose coordinates are specified by these two flags. The default is the origin `0 + 0i`.
    * `-j` enable Julia mode  
    * `-A LNG_DBL` sets real part of the constant c in Julia mode
    * `-B LNG_DBL` sets imaginary part of the constant c in Julia mode  
    	In Julia mode the constant `c` in the iterative formulas to render the different fractal is fixed, and doesn't depend on the starting point in the complex plane. These three flags allow the user to run the rendering in Julia 	mode and to set the real and imaginary part of this constant. The default `c` is the origin `0 + 0i`.
    * `-s LNG_DBL` sets scaling factor. The default is 1.  
    	Changes the amount of zoom done around the centeral rendering point.
    * `-t INT` sets maximum number of iterations. The default is 2000.
    * `-b LNG_DBL` sets bailout radius for the fractal. The default is 2.
    * `-T INT` sets maximum number of threads working on the rendering of the image. The default is 8.
* **Coloring related flags**
	* `-p STRING` sets filename of input config file for color palette to the content of `STRING`.  
		When this flag is used, the program looks in the folder into which the executable is placed, for
        a file whose name is the content of `STRING`. It must be a simple plain text file.  
        Each line in this file represents a color that is then loaded in the color palette. The colors must be specified by writing their individual `R`(ed), `G`(reen) and `B`(lue) channel values separated by a space.  
        For example, a line containing `255 255 255` represents white, `0 0 0` represents black and `255 255 0` is
        yellow.  
        In general, the format for a single line is `uint8_t uint8_t uint8_t`.  
        By default, the program looks for a file named `palette` (note that the name has no extension. If you want to load `palette.txt` you will have to specify it) and if the palette file is not found or an error occurs while trying to load the colors, the palette is loaded with only black and white.
    * `-c INT` sets coloring mode. Must be a number in [0, 4]. The default is 3, because it's the best looking one.  
    	The supported coloring modes are:
    	* `0` - binary (black or white)  
        * `1` - ln (which is the most popular option but I'm really struggling in implementing it properly. It doesn't look great)  
        * `2` - square root of iterations  
        * `3` - cubic root of iterations  
        * `4` - angle of the last iterated point (calculated through `atan2`)
* **Other flags**
	* `-v`print verbose output, more user friendly.  
		This flag enables the printing to `stdout` of some useful info about the parameters that have been set and the state of the rendering.
        
## Automation
This program has been thought to be easily integrated with other programs.  
An example is the `bash` script called `zoom.sh`, which repeatedly runs `liten-fraktaler` at different levels of zoom to generate a zoom in sequence to a Minibrot in the [Mandelbrot set](https://en.wikipedia.org/wiki/Mandelbrot_set) located at `-1.992561950824231 + 0i`. This script requires [bc](https://www.gnu.org/software/bc/) to work; it can be easily modified to work with other fractals and do other types of zoom/pan sequences.

## To do
A lot of stuff, there are some sparse `TODO` in the source code, but summarized here they are:
* add multi letter and multi parameter flags
* add flag to change the fractal rotation in the plane
* plane warping
* add flag to automatically set the number of iterations based on the specified zoom
* [Pickover stalks](https://en.wikipedia.org/wiki/Pickover_stalk)
* allow image segmentation, i.e. subdivide big images that wouldn't fit in RAM, into multiple smaller images, compute each of the individually, and save them on disk
* fix the logarithmic coloring mode and add many others
* fix the spade fractal created by [ahhhh6980](https://github.com/ahhhh6980), the output image of my program is very different from his and I don't know why
* add more fractals
* reorganize stuff
* optimize stuff
* arbitrary precision arithmetics

So yeah, if you want to contribute you're more than welcome! ;)

## Working on:
* adding fractals
* reorganize stuff

Update with lots of new fractals coming soon!
