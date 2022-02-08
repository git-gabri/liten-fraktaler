# Liten fraktaler
A fractal rendering software, written in C++, inspired by a much more popular program called [Kalles Fraktaler](https://mathr.co.uk/kf/kf.html)

## Compiling on Linux with `cmake` and `gcc`
* Clone this repository
* `cd` into the repository where the `CMakeLists.txt` is located
* `cmake -DCMAKE_EXPORT_COMPILE_COMMANDS:BOOL=TRUE -DCMAKE_BUILD_TYPE:STRING=Release -B./build`
* `cmake --build ./build/ --config Release --target all`

You should now have an executable in the `build` folder called `liten-fraktaler`.

## Usage
The help is available by launching the program with the `-H` or the `--help` flag, like this:  
`./liten-fraktaler --help`.
        
## Automation
This program has been thought to be easily integrated with other programs.  
An example is the `bash` script called `zoom.sh`, which repeatedly runs `liten_fraktaler` at different levels of zoom to generate a zoom in sequence to a Minibrot in the [Mandelbrot set](https://en.wikipedia.org/wiki/Mandelbrot_set) located at `-1.992561950824231 + 0i`. This script requires [bc](https://www.gnu.org/software/bc/) to work; it can be easily modified to work with other fractals and do other types of zoom/pan sequences.

## To do
Still a lot of stuff, there are some sparse `TODOs` in the source code (and some personal notes in the `notes` file), but the doable ones are summarized here:
* fix smooth coloring
* solid guessing
* supersampling
* add flag to automatically set the number of iterations based on the specified zoom
* [Pickover stalks](https://en.wikipedia.org/wiki/Pickover_stalk)
* allow image segmentation, i.e. subdivide big images that wouldn't fit in RAM, into multiple smaller images, compute each of the individually, and save them on disk
* hybrid fractals
* reorganize and refactor
* [...]